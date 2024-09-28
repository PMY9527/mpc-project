#include "FSM/Mpc.h"
#include <iomanip>

Mpc::Mpc(CtrlComponents *ctrlComp)
             :FSMState(ctrlComp, FSMStateName::TROTTING, "trottingMPC"), 
              _est(ctrlComp->estimator), _phase(ctrlComp->phase), 
              _contact(ctrlComp->contact), _robModel(ctrlComp->robotModel), 
              _balCtrl(ctrlComp->balCtrl){
    _gait = new GaitGenerator(ctrlComp); // GaitGenerator类的指针，生成期望位置和速度：_pos-FeetGlobalGoal 和 _velFeetGlobalGoal。

    _gaitHeight = 0.08;

    Mpc::~Mpc() {
    }

#ifdef ROBOT_TYPE_Go1
    _Kpp = Vec3(70, 70, 70).asDiagonal();
    _Kdp = Vec3(10, 10, 10).asDiagonal();
    _kpw = 780; 
    _Kdw = Vec3(70, 70, 70).asDiagonal();
    _KpSwing = Vec3(400, 400, 400).asDiagonal();
    _KdSwing = Vec3(10, 10, 10).asDiagonal();
#endif

#ifdef ROBOT_TYPE_A1
    _Kpp = Vec3(20, 20, 100).asDiagonal();
    _Kdp = Vec3(20, 20, 20).asDiagonal();
    _kpw = 400;
    _Kdw = Vec3(50, 50, 50).asDiagonal();
    _KpSwing = Vec3(400, 400, 400).asDiagonal();
    _KdSwing = Vec3(10, 10, 10).asDiagonal();
#endif

    _vxLim = _robModel->getRobVelLimitX();
    _vyLim = _robModel->getRobVelLimitY();
    _wyawLim = _robModel->getRobVelLimitYaw();

}

Mpc::~Mpc(){
    delete _gait;
}

void Mpc::enter(){ // 初始化
    
    _pcd = _est->getPosition(); 
    _pcd(2) = -_robModel->getFeetPosIdeal()(2, 0); 
    _vCmdBody.setZero();
    _yawCmd = _lowState->getYaw();
    _Rd = rotz(_yawCmd);
    _wCmdGlobal.setZero();

    _ctrlComp->ioInter->zeroCmdPanel();
    _gait->restart();
}

void Mpc::exit(){ 
    _ctrlComp->ioInter->zeroCmdPanel();
    _ctrlComp->setAllSwing(); 
}

FSMStateName Mpc::checkChange(){
    if(_lowState->userCmd == UserCommand::L2_B){
        return FSMStateName::PASSIVE;
    }
    else if(_lowState->userCmd == UserCommand::L2_A){
        return FSMStateName::FIXEDSTAND;
    }
    else{
        return FSMStateName::TROTTING;
    }
}

void Mpc::run(){ // here
    _posBody = _est->getPosition(); // 现时的估计状态值（位置）， Eigen::Matrix<double, 3, 1>
    _velBody = _est->getVelocity(); // 现时的估计状态值（速度）， Eigen::Matrix<double, 3, 1>
    _posFeet2BGlobal = _est->getPosFeet2BGlobal(); // 足端在世界坐标系下相对于机身中心的位置坐标。Eigen::Matrix<double, 3, 4>
    _posFeetGlobal = _est->getFeetPos(); // 足端在s下的位置坐标。（四只脚--3by4）Eigen::Matrix<double, 3, 4>
    _velFeetGlobal = _est->getFeetVel(); // 足端在s下的速度。Eigen::Matrix<double, 3, 4>
    _mass = robModel->getRobMass(); // 重量 
    R_sb = _lowState->getRotMat(); // R_sb 机身姿态，世界 -> 机身。MAT3
    R_bs = R_sb.transpose(); // R_bs , 机身 -> 世界。MAT3
    _yaw = _lowState->getYaw(); // 机身偏航角，about z axis
    _dYaw = _lowState->getDYaw(); // 机身偏航速度
    _userValue = _lowState->userValue; 

    getUserCmd(); 
    calcCmd();

    getAB();

    _gait->setGait(_vCmdGlobal.segment(0,2), _wCmdGlobal(2), _gaitHeight);
    _gait->run(_posFeetGlobalGoal, _velFeetGlobalGoal);

    for(int i(0); i<4; ++i){
        if((*_contact)(i) == 0){
            _lowCmd->setSwingGain(i);
        }else{
            _lowCmd->setStableGain(i);
        }
    }

}

bool Mpc::checkStepOrNot(){
    if( (fabs(_vCmdBody(0)) > 0.03) ||
        (fabs(_vCmdBody(1)) > 0.03) ||
        (fabs(_posError(0)) > 0.08) ||
        (fabs(_posError(1)) > 0.08) ||
        (fabs(_velError(0)) > 0.05) ||
        (fabs(_velError(1)) > 0.05) ||
        (fabs(_dYawCmd) > 0.20) ){
        return true;
    }else{
        return false;
    }
}

void Mpc::setHighCmd(double vx, double vy, double wz){ //here
    _vCmdBody(0) = vx;
    _vCmdBody(1) = vy;
    _vCmdBody(2) = 0; 
    _dYawCmd = wz;
}

void Mpc::getUserCmd(){ // 获取用户指令。
    /* Movement */
    _vCmdBody(0) =  invNormalize(_userValue.ly, _vxLim(0), _vxLim(1));
    _vCmdBody(1) = -invNormalize(_userValue.lx, _vyLim(0), _vyLim(1));
    _vCmdBody(2) = 0; 

    /* Turning */
    _dYawCmd = -invNormalize(_userValue.rx, _wyawLim(0), _wyawLim(1));
    _dYawCmd = 0.9*_dYawCmdPast + (1-0.9) * _dYawCmd;
    _dYawCmdPast = _dYawCmd;
}

void Mpc::calcCmd(){ // 将用户指令换算成期望线速度和角速度。
    /* Movement */
    _vCmdGlobal = R_sb * _vCmdBody;

    _vCmdGlobal(0) = saturation(_vCmdGlobal(0), Vec2(_velBody(0)-0.2, _velBody(0)+0.2));
    _vCmdGlobal(1) = saturation(_vCmdGlobal(1), Vec2(_velBody(1)-0.2, _velBody(1)+0.2));

    _pcd(0) = saturation(_pcd(0) + _vCmdGlobal(0) * _ctrlComp->dt, Vec2(_posBody(0) - 0.05, _posBody(0) + 0.05));
    _pcd(1) = saturation(_pcd(1) + _vCmdGlobal(1) * _ctrlComp->dt, Vec2(_posBody(1) - 0.05, _posBody(1) + 0.05));

    _vCmdGlobal(2) = 0;

    /* Turning */
    _yawCmd = _yawCmd + _dYawCmd * _ctrlComp->dt;

    _Rd = rotz(_yawCmd);
    _wCmdGlobal(2) = _dYawCmd;
}

void Mpc::getAB(){ //here
    
    I_Body(0,0) = 0.1320;
    I_Body(1,1) = 0.3475;
    I_Body(2,2) = 0.3775; // 机身坐标下，a1的惯性张量。
    I_World = R_sb * I_Body * R_bs;
    I_WorldINV = I_World.inverse();

    _yaw = _lowState->getYaw();

    R_yaw << std::cos(_yaw), -std::sin(_yaw), 0,
             std::sin(_yaw),  std::cos(_yaw), 0,
             0,               0,              1;

    _continuousA.block(0,6,3,3) = R_yaw.transpose();
    _continuousA.block(3,9,3,3) = I3;
    _continuousA.block(9,12,3,1) = Eigen::Vector3d(0, 0, 1);
    _continuousA.block(12,12,1,1) = 1;
    _discretizedA = I13 + ctrlComp->dt * _continuousA;

    _posFeet2BGlobal = _est->getPosFeet2BGlobal();

    for(int i = 0; i < 4; i++){
        _discretizedB.block(6,i*3,3,3) = ctrlComp->dt * skew( _posFeet2BGlobal.col(i)) * I_WorldINV;  
        _discretizedB.block(9,i*3,3,3) = I3 * ctrlComp->dt / _mass;
    }
}

void Mpc::SetupCostfunc(){
    int N = 10; // Prediction Horizon
    std::vector<Mat> A_powers(N);

    A_powers[0] = _discretizedA;
    for (int i = 1; i < N; ++i) {
        A_powers[i] = A_powers[i - 1] * (_discretizedA);
    }
    
    Mat Aqp = Mat::Zero(N * 13, 13);
    for (int i = 0; i < N; ++i) {
        Aqp.block(i * 13, 0, 13, 13) = A_powers[i];
    }

    Mat Bqp = Mat::Zero(N * 13, N * 12);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            Mat A_power = (i - j - 1 >= 0) ? A_powers[i - j - 1] : I13;
            Bqp.block(i * 13, j * 12, 13, 12) = A_power * _discretizedB;
        }
    }

    Mat Qq = Mat::Identity(N * 13, N * 13); // 状态权重
    Mat Rr = Mat::Identity(N * 12, N * 12); // 控制输入权重
    Mat Q = 5 * Qq;
    Mat R = 1 * Rr;
    
    Vec _CurrentStatesVec(13);
    _CurrentStatesVec <<               0,           0,        _yaw,
                             _posBody(0), _posBody(1), _posBody(2),
                                       0,           0,       _dYaw,
                             _velBody(0), _velBody(1), _velBody(2),
                                   -9.81;
    
    Vec X_ref = Vec::Zero(N * 13);

    Vec _TargetStatesVec(13);
    _TargetStatesVec <<          0,              0,        _yawCmd,
                           _pcd(0),        _pcd(1),        _pcd(2),
                                 0,              0,       _dYawCmd,
                    _vCmdGlobal(0), _vCmdGlobal(1), _vCmdGlobal(2),
                             -9.81;

    for (int i = 0; i < N; ++i){ 
     X_ref.segment(i * 13, 13) = _TargetStatesVec;
    }

    Mat Hqp = 2 * Bqp.transpose() * Q * Bqp + R;
    Mat fqp = 2 * Bqp.transpose() * Q * (Aqp *_CurrentStatesVec - X_ref);
    
} 

void Mpc::ConstraintSetup(){
    // Constraints (friction cone) ：
            // 摆动腿的解（u / 力） = 0。
            // 触地腿不能打滑，即u（地面反作用力）的水平分量不能超过f_max。

    double mu = 0.4;
    
    Cons_singleLEG <<  -1,  0, mu, // Single Leg Constraints
                       0,  -1, mu,
                       1,  mu, mu,
                       0,   1, mu,
                       0,   0,  1;
    
    
    for(int i = 0; i < N * 4; ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
     fmat.block(i*5,i*3,5,3) = Cons_singleLEG; // 触地腿的约束。
   }
     
    for(int i(0); i<4; ++i){
        if((*_contact)(j) == 0){  

        }
    }


}

void Mpc::solveQP(){}



