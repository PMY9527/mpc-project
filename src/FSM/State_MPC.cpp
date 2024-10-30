/**********************************************************************
 Copyright (c) 2020-2023, Unitree Robotics.Co.Ltd. All rights reserved.
***********************************************************************/
#include "FSM/State_MPC.h"
#include <iomanip>

State_MPC::State_MPC(CtrlComponents *ctrlComp)
    : FSMState(ctrlComp, FSMStateName::MPC, "mpc"),
      _est(ctrlComp->estimator), _phase(ctrlComp->phase),
      _contact(ctrlComp->contact), _robModel(ctrlComp->robotModel),
      _balCtrl(ctrlComp->balCtrl)
{
    t1_prev = std::chrono::high_resolution_clock::now(); // initialise the timer for MPC benchmarking
    _gait = new GaitGenerator(ctrlComp); 
    _gaitHeight = 0.08;

    // this is for unitree A1
    _Kpp = Vec3(20, 20, 100).asDiagonal(); 
    _Kdp = Vec3(20, 20, 20).asDiagonal();
    _kpw = 400;
    _Kdw = Vec3(50, 50, 50).asDiagonal();
    _KpSwing = Vec3(400, 400, 400).asDiagonal();
    _KdSwing = Vec3(10, 10, 10).asDiagonal();

    _vxLim = _robModel->getRobVelLimitX();
    _vyLim = _robModel->getRobVelLimitY();
    _wyawLim = _robModel->getRobVelLimitYaw();

    _mass = _robModel->getRobMass();
    Ic = _robModel->getRobInertial();

    max[0] = 1000.0;
    max[1] = 1000.0;
    max[2] = -3.0 * _mass * g;

    min[0] = -1000.0;
    min[1] = -1000.0;
    min[2] = 0.0;

    miuMat << 1, 0, miu,
        -1, 0, miu,
        0, 1, miu,
        0, -1, miu,
        0, 0, 1;

    setWeight();
}

void State_MPC::setWeight()
{
    Eigen::Matrix<double, 1, nx> Q_diag;
    Eigen::Matrix<double, 1, nu> R_diag;

    Q_diag << 270.0, 270.0, 1.0, // 欧拉角(roll pitch yaw)
            5.0, 5.0, 270.0, // pCoM
            1.0, 1.0, 20.0, // w
            20.0, 20.0, 40.0, //vcom 
        
    R_diag <<   1.0, 1.0, 0.1, // 右前 (xyz)
                1.0, 1.0, 0.1, // 左前
                1.0, 1.0, 0.1, // 右后
                1.0, 1.0, 0.1, // 左后
                0.0;
    R_diag = R_diag * 0.8 * 1e-3; 

    Eigen::MatrixXd Q_diag_N = Eigen::MatrixXd::Zero(1, nx * mpc_N);
    Eigen::MatrixXd R_diag_N = Eigen::MatrixXd::Zero(1, nu * ch);

    Q = Eigen::MatrixXd::Zero(nx * mpc_N, nx * mpc_N);
    R = Eigen::MatrixXd::Zero(nu * ch, nu * ch);

    for (int i = 0; i < mpc_N; i++)
    {
        Q_diag_N.block<1, nx>(0, i * nx) = Q_diag;
    }

    for (int i = 0; i < ch; i++)
    {
        R_diag_N.block<1, nu>(0, i * nu) = R_diag;
    }

    for (int i = 0; i < nx * mpc_N; i++)
    {
        Q(i, i) = Q_diag_N(0, i);
    }

    for (int i = 0; i < nu * ch; i++)
    {
        R(i, i) = R_diag_N(0, i);
    }
}

State_MPC::~State_MPC()
{
    delete _gait;
}

void State_MPC::enter()
{
    _pcd = _est->getPosition();
    _pcd(2) = -_robModel->getFeetPosIdeal()(2, 0);
    _vCmdBody.setZero();
    _yawCmd = _lowState->getYaw();
    _Rd = rotz(_yawCmd);
    _wCmdGlobal.setZero();
    _ctrlComp->ioInter->zeroCmdPanel();
    _gait->restart();

    _forceFeetGlobal.setZero();
    _tau.setZero();
}

void State_MPC::exit()
{
    _ctrlComp->ioInter->zeroCmdPanel();
    _ctrlComp->setAllSwing();
}

FSMStateName State_MPC::checkChange()
{
    if (_lowState->userCmd == UserCommand::L2_B || (_forceFeetGlobal.array() != _forceFeetGlobal.array()).any()) // if nan, quit to passive
    {
        return FSMStateName::PASSIVE;
    }
    else if (_lowState->userCmd == UserCommand::L2_A)
    {
        return FSMStateName::FIXEDSTAND;
    }
    else
    {
        return FSMStateName::MPC;
    }
}

void State_MPC::run()
{   
    _posBody = _est->getPosition();
    _velBody = _est->getVelocity();
    _posFeet2BGlobal = _est->getPosFeet2BGlobal();
    _posFeetGlobal = _est->getFeetPos();
    _velFeetGlobal = _est->getFeetVel();
    _B2G_RotMat = _lowState->getRotMat();
    _G2B_RotMat = _B2G_RotMat.transpose();
    _yaw = _lowState->getYaw();
    _dYaw = _lowState->getDYaw();

    _userValue = _lowState->userValue;

    getUserCmd();
    calcCmd();

    _gait->setGait(_vCmdGlobal.segment(0, 2), _wCmdGlobal(2), _gaitHeight);
    _gait->run(_posFeetGlobalGoal, _velFeetGlobalGoal);

    calcTau();
    calcQQd(); // q and qd 

    _ctrlComp->setStartWave();

    _lowCmd->setTau(_tau); // tau limit +_50 
    _lowCmd->setQ(vec34ToVec12(_qGoal));
    _lowCmd->setQd(vec34ToVec12(_qdGoal));

    for (int i(0); i < 4; ++i)
    {
        if ((*_contact)(i) == 0)
        {
            _lowCmd->setSwingGain(i);
        }
        else
        {
            _lowCmd->setStableGain(i);
        }
    }
}

void State_MPC::setHighCmd(double vx, double vy, double wz)
{
    _vCmdBody(0) = vx;
    _vCmdBody(1) = vy;
    _vCmdBody(2) = 0;
    _dYawCmd = wz;
}

void State_MPC::getUserCmd()
{
    /* Movement */
    _vCmdBody(0) = invNormalize(_userValue.ly, _vxLim(0), _vxLim(1));
    _vCmdBody(1) = -invNormalize(_userValue.lx, _vyLim(0), _vyLim(1));
    _vCmdBody(2) = 0;

    /* Turning */
    _dYawCmd = -invNormalize(_userValue.rx, _wyawLim(0), _wyawLim(1));
    _dYawCmd = 0.9 * _dYawCmdPast + (1 - 0.9) * _dYawCmd;
    _dYawCmdPast = _dYawCmd;
}

void State_MPC::calcCmd()
{
    /* Movement */
    _vCmdGlobal = _B2G_RotMat * _vCmdBody;

    _vCmdGlobal(0) = saturation(_vCmdGlobal(0), Vec2(_velBody(0) - 0.2, _velBody(0) + 0.2));
    _vCmdGlobal(1) = saturation(_vCmdGlobal(1), Vec2(_velBody(1) - 0.2, _velBody(1) + 0.2));

    _pcd(0) = saturation(_pcd(0) + _vCmdGlobal(0) * _ctrlComp->dt, Vec2(_posBody(0) - 0.05, _posBody(0) + 0.05));
    _pcd(1) = saturation(_pcd(1) + _vCmdGlobal(1) * _ctrlComp->dt, Vec2(_posBody(1) - 0.05, _posBody(1) + 0.05));

    _vCmdGlobal(2) = 0;

    /* Turning */
    _yawCmd = _yawCmd + _dYawCmd * _ctrlComp->dt;

    _Rd = rotz(_yawCmd);
    _wCmdGlobal(2) = _dYawCmd;
}

void State_MPC::calcTau()
{
    _posError = _pcd - _posBody;
    _velError = _vCmdGlobal - _velBody;

    _ddPcd = _Kpp * _posError + _Kdp * _velError;
    _dWbd = _kpw * rotMatToExp(_Rd * _G2B_RotMat) + _Kdw * (_wCmdGlobal - _lowState->getGyroGlobal());

    _ddPcd(0) = saturation(_ddPcd(0), Vec2(-3, 3));
    _ddPcd(1) = saturation(_ddPcd(1), Vec2(-3, 3));
    _ddPcd(2) = saturation(_ddPcd(2), Vec2(-5, 5));

    _dWbd(0) = saturation(_dWbd(0), Vec2(-40, 40));
    _dWbd(1) = saturation(_dWbd(1), Vec2(-40, 40));
    _dWbd(2) = saturation(_dWbd(2), Vec2(-10, 10));

    calcFe();

    std::cout << "********forceFeetGlobal(MPC)********" << std::endl
              << _forceFeetGlobal << std::endl;
            
    //std::cout << "********Tau(MPC)********" << std::endl
    //          << _tau << std::endl;


    for (int i(0); i < 4; ++i)
    {
        if ((*_contact)(i) == 0) // 摆动腿
        { 
            _forceFeetGlobal.col(i) = _KpSwing * (_posFeetGlobalGoal.col(i) - _posFeetGlobal.col(i)) + _KdSwing * (_velFeetGlobalGoal.col(i) - _velFeetGlobal.col(i));
        }
    }

    //std::cout << "********_contact(MPC)********" << std::endl;
    //for (int i = 0; i < 4; ++i) {
    //    std::cout << (*_contact)(i) << " ";
    //}
    //std::cout << std::endl;


    _forceFeetBody = _G2B_RotMat * _forceFeetGlobal;
    _q = vec34ToVec12(_lowState->getQ());
    _tau = _robModel->getTau(_q, _forceFeetBody);

}

void State_MPC::calcQQd()
{

    Vec34 _posFeet2B;
    _posFeet2B = _robModel->getFeet2BPositions(*_lowState, FrameType::BODY);

    for (int i(0); i < 4; ++i)
    {
        _posFeet2BGoal.col(i) = _G2B_RotMat * (_posFeetGlobalGoal.col(i) - _posBody);
        _velFeet2BGoal.col(i) = _G2B_RotMat * (_velFeetGlobalGoal.col(i) - _velBody);
        // _velFeet2BGoal.col(i) = _G2B_RotMat * (_velFeetGlobalGoal.col(i) - _velBody - _B2G_RotMat * (skew(_lowState->getGyro()) * _posFeet2B.col(i)) );  //  c.f formula (6.12)
    }

    _qGoal = vec12ToVec34(_robModel->getQ(_posFeet2BGoal, FrameType::BODY)); // _qGoal 各个关节的目标角度
    _qdGoal = vec12ToVec34(_robModel->getQd(_posFeet2B, _velFeet2BGoal, FrameType::BODY)); // _qdGoal 各个关节的目标角速度
}

#undef inverse
void State_MPC::calcFe()
{
    // 当前状态：欧拉角、机身位置、角速度、机身速度
    X_cur << _G2B_RotMat.eulerAngles(2, 1, 0), _posBody, _lowState->getGyroGlobal(), _velBody;
   // std::cout << "X_cur" << std::endl 
   //          << X_cur << std::endl;

    // 设置期望状态
    for (int i = 0; i < (mpc_N - 1); i++)
        Xd.block<nx, 1>(nx * i, 0) = Xd.block<nx, 1>(nx * (i + 1), 0);
        Xd(nx * (mpc_N - 1) + 2) = _yawCmd;
    for (int j = 0; j < 3; j++)
        Xd(nx * (mpc_N - 1) + 3 + j) = _pcd(j);
    for (int j = 0; j < 3; j++)
        Xd(nx * (mpc_N - 1) + 6 + j) = _wCmdGlobal(j);
    for (int j = 0; j < 3; j++)
        Xd(nx * (mpc_N - 1) + 9 + j) = _vCmdGlobal(j);

    // 构造Ac

    for (int i = 0; i < mpc_N; i++)
    {
        R_curz[i] = Rz3(_yaw);
        Ac[i].setZero();
        Ac[i].block<3, 3>(0, 6) = R_curz[i].transpose();
        Ac[i].block<3, 3>(3, 9) = Eigen::MatrixXd::Identity(3, 3);
        A[i] = Eigen::MatrixXd::Identity(nx, nx) + _ctrlComp->dt * Ac[i];
    }

    // 构造Bc
    for (int i = 0; i < mpc_N; i++)
    {
        Mat3 Ic_W_inv;
        Ic_W_inv = (R_curz[i] * Ic * R_curz[i].transpose()).inverse();
        Bc[i].setZero();
        Bc[i].block<3, 3>(6, 0) = Ic_W_inv * CrossProduct_A(_posFeet2BGlobal.block<3, 1>(0, 0));
        Bc[i].block<3, 3>(6, 3) = Ic_W_inv * CrossProduct_A(_posFeet2BGlobal.block<3, 1>(0, 1));
        Bc[i].block<3, 3>(6, 6) = Ic_W_inv * CrossProduct_A(_posFeet2BGlobal.block<3, 1>(0, 2));
        Bc[i].block<3, 3>(6, 9) = Ic_W_inv * CrossProduct_A(_posFeet2BGlobal.block<3, 1>(0, 3));
        Bc[i].block<3, 3>(9, 0) = I3 / _mass;
        Bc[i].block<3, 3>(9, 3) = I3 / _mass;
        Bc[i].block<3, 3>(9, 6) = I3 / _mass;
        Bc[i].block<3, 3>(9, 9) = I3 / _mass;
        Bc[i]((nx - 1), (nu - 1)) = 1.0;
        B[i] = _ctrlComp->dt * Bc[i];
    }

    // 构造Aqp
    Aqp.setZero();
    for (int i = 0; i < mpc_N; i++)
        Aqp.block<nx, nx>(i * nx, 0) = Eigen::MatrixXd::Identity(nx, nx);
    for (int i = 0; i < mpc_N; i++)
        for (int j = 0; j < i + 1; j++)
            Aqp.block<nx, nx>(i * nx, 0) = A[j] * Aqp.block<nx, nx>(i * nx, 0);

    // 构造Bqp
    Aqp1.setZero();
    Bqp1.setZero();
    Bqp.setZero();
    for (int i = 0; i < mpc_N; i++)
        for (int j = 0; j < i + 1; j++)
            Aqp1.block<nx, nx>(i * nx, j * nx) = Eigen::MatrixXd::Identity(nx, nx);
    for (int i = 1; i < mpc_N; i++)
        for (int j = 0; j < i; j++)
            for (int k = j + 1; k < (i + 1); k++)
                Aqp1.block<nx, nx>(i * nx, j * nx) = A[k] * Aqp1.block<nx, nx>(i * nx, j * nx);

    for (int i = 0; i < mpc_N; i++)
        Bqp1.block<nx, nu>(i * nx, i * nu) = B[i];
    Eigen::MatrixXd Bqp11 = Eigen::MatrixXd::Zero(nu * mpc_N, nu * ch);
    Bqp11.setZero();
    Bqp11.block<nu * ch, nu * ch>(0, 0) = Eigen::MatrixXd::Identity(nu * ch, nu * ch);
    for (int i = 0; i < (mpc_N - ch); i++)
        Bqp11.block<nu, nu>(nu * ch + i * nu, nu * (ch - 1)) = Eigen::MatrixXd::Identity(nu, nu);

    Eigen::MatrixXd B_tmp = Eigen::MatrixXd::Zero(nx * mpc_N, nu * ch);
    B_tmp = Bqp1 * Bqp11;
    Bqp = Aqp1 * B_tmp;

    // 二次规划一般形式 1/2 x^T*H*x + x^T * c
    H_ = 2 * (Bqp.transpose() * Q * Bqp + R);
    c_ = 2 * Bqp.transpose() * Q * (Aqp * X_cur - Xd);

    Constraints();
    solveQP();

    _forceFeetGlobal = vec12ToVec34(F_);
}

void State_MPC::Constraints()
{
    int contactLegNum = 0;
    for (int i(0); i < 4; ++i)
    {
        if ((*_contact)(i) == 1)
        {
            contactLegNum += 1;
        }
    }

    CI_.resize(5 * contactLegNum * ch, nu * ch);
    ci0_.resize(5 * contactLegNum * ch);
    CE_.resize((3 * (4 - contactLegNum) + 1) * ch, nu * ch);
    ce0_.resize((3 * (4 - contactLegNum) + 1) * ch);

    CI_.setZero();
    ci0_.setZero();
    CE_.setZero();
    ce0_.setZero();

    for (int k = 0; k < ch; k++)
    {
        int ceID = 0;
        int ciID = 0;
        for (int i(0); i < 4; ++i)
        {
            if ((*_contact)(i) == 1)
            {
                CI_.block<5, 3>(5 * contactLegNum * k + 5 * ciID, nu * k + 3 * i) = miuMat;
                ++ciID;
            }
            else
            {
                CE_.block<3, 3>((3 * (4 - contactLegNum) + 1) * k + 3 * ceID, nu * k + 3 * i) = I3;
                ++ceID;
            }
        }
        CE_((3 * (4 - contactLegNum) + 1) * k + 3 * (4 - contactLegNum), nu * k + nu - 1) = 1.0;
        ce0_[(3 * (4 - contactLegNum) + 1) * k + 3 * (4 - contactLegNum)] = -g;
    }
}

void State_MPC::solveQP()
{
    int n = nu * ch;
    int m = ce0_.size();
    int p = ci0_.size();

    G.resize(n, n);
    CE.resize(n, m);
    CI.resize(n, p);
    g0.resize(n);
    ce0.resize(m);
    ci0.resize(p);
    x.resize(n);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            G[i][j] = H_(i, j);
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            CE[i][j] = (CE_.transpose())(i, j);
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            CI[i][j] = (CI_.transpose())(i, j);
        }
    }

    for (int i = 0; i < n; ++i)
    {
        g0[i] = c_[i];
    }

    for (int i = 0; i < m; ++i)
    {
        ce0[i] = ce0_[i];
    }

    for (int i = 0; i < p; ++i)
    {
        ci0[i] = ci0_[i];
    }

    //std::cout << "n:" << n << "m:" << m << "p:" << p << std::endl;
    double value = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double_1 = t1 - t1_prev;
    t1_prev = t1;

    std::cout << "time for each mpc loop: " << ms_double_1.count() << "ms" << std::endl; // around 2.5 ms
    
    // it runs for about 0.1 to 0.2 ms, meaning 
    
    for (int i = 0; i < 12; ++i)
    {
        F_[i] = -x[i];
    }
}

Eigen::Matrix<double, 3, 3> State_MPC::CrossProduct_A(Eigen::Matrix<double, 3, 1> A)
{
    Eigen::Matrix<double, 3, 3> M;
    M << 0.0, -A[2], A[1],
        A[2], 0.0, -A[0],
        -A[1], A[0], 0.0;
    return M;
}

Eigen::Matrix<double, 3, 3> State_MPC::Rz3(double theta)
{ // local to world
    // for 2D-XY vector, rotation matrix along z axis
    Eigen::Matrix<double, 3, 3> M;
    M << cos(theta), -sin(theta), 0,
        sin(theta), cos(theta), 0,
        0, 0, 1;
    return M;
}