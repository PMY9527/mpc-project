/**********************************************************************
 Copyright (c) 2020-2023, Unitree Robotics.Co.Ltd. All rights reserved.
***********************************************************************/
#ifndef MPC_H
#define MPC_H

#include "FSM/FSMState.h"
#include "Gait/GaitGenerator.h"
#include "control/BalanceCtrl.h"
#include "thirdParty/quadProgpp/QuadProg++.hh"
#include "thirdParty/quadProgpp/Array.hh"
#include <chrono>
#include <vector>
#include <iostream>


static const int mpc_N = 10; // 计算mpc_N步状态
static const int ch = 3;     // 为了让Bqp的维数不那么庞大，减小运算量，将Bqp的列数由13*10减小为13*3
static const int nx = 12;    // 状态向量的维数
static const int nu = 13;    // 控制输入的维数
static const int nc = 20;    // 约束的维数

static const double NEGATIVE_NUMBER = -1000000.0; // 负无穷
static const double POSITIVE_NUMBER = 1000000.0;  // 正无穷

static const double g = -9.8;
static const double miu = 0.5; // 摩擦系数

class State_MPC : public FSMState
{
public:
    State_MPC(CtrlComponents *ctrlComp);
    ~State_MPC();
    void enter();
    void run();
    void exit();
    virtual FSMStateName checkChange();
    void setHighCmd(double vx, double vy, double wz);

private:
    void calcTau();
    void calcQQd();
    void calcCmd();
    virtual void getUserCmd();
    bool checkStepOrNot();
    void calcFe();

    GaitGenerator *_gait;
    Estimator *_est;
    QuadrupedRobot *_robModel;
    BalanceCtrl *_balCtrl;

    // Rob State
    Vec3 _posBody, _velBody;
    double _yaw, _dYaw;
    Vec34 _posFeetGlobal, _velFeetGlobal;
    Vec34 _posFeet2BGlobal;
    RotMat _B2G_RotMat, _G2B_RotMat;
    Vec12 _q;

    // Robot command
    Vec3 _pcd;
    Vec3 _vCmdGlobal, _vCmdBody;
    double _yawCmd, _dYawCmd;
    double _dYawCmdPast;
    Vec3 _wCmdGlobal;
    Vec34 _posFeetGlobalGoal, _velFeetGlobalGoal;
    Vec34 _posFeet2BGoal, _velFeet2BGoal;
    RotMat _Rd;
    Vec3 _ddPcd, _dWbd;
    Vec34 _forceFeetGlobal, _forceFeetBody;
    Vec34 _qGoal, _qdGoal;
    Vec12 _tau, tau_set;

    // Control Parameters
    double _gaitHeight;
    Vec3 _posError, _velError;
    Mat3 _Kpp, _Kdp, _Kdw;
    double _kpw;
    Mat3 _KpSwing, _KdSwing;
    Vec2 _vxLim, _vyLim, _wyawLim;
    Vec4 *_phase;
    VecInt4 *_contact;

    // Calculate average value
    AvgCov *_avg_posError = new AvgCov(3, "_posError", true, 1000, 1000, 1);
    AvgCov *_avg_angError = new AvgCov(3, "_angError", true, 1000, 1000, 1000);

    // MPC用到的变量
    double _mass;          // 质量
    double max[3], min[3]; // 足底力最大值最小值
    Eigen::Matrix<double, 5, 3> miuMat;

    Eigen::Matrix<double, 3, 3> Ic;
    Eigen::Matrix<double, nx, 1> X_cur;
    Eigen::Matrix<double, nx * mpc_N, 1> Xd;
    Eigen::Matrix<double, 3, 3> R_curz[mpc_N];
    Eigen::Matrix<double, nx, nx> Ac[mpc_N], A[mpc_N];
    Eigen::Matrix<double, nx, nu> Bc[mpc_N], B[mpc_N];

    Eigen::Matrix<double, nx * mpc_N, nx> Aqp;
    Eigen::Matrix<double, nx * mpc_N, nx * mpc_N> Aqp1;
    Eigen::Matrix<double, nx * mpc_N, nu * mpc_N> Bqp1;
    Eigen::Matrix<double, nx * mpc_N, nu * ch> Bqp;
    Eigen::Matrix<double, nu * ch, nu * ch> H_;
    Eigen::Matrix<double, nu * ch, 1> c_;
    Eigen::Matrix<double, nx * mpc_N, nx * mpc_N> Q;
    Eigen::Matrix<double, nu * ch, nu * ch> R;
    Eigen::Matrix<double, nu * ch, 1> lb_qp;
    Eigen::Matrix<double, nu * ch, 1> ub_qp;
    Eigen::Matrix<double, nc, nu> C_control;
    Eigen::Matrix<double, nc * ch, nu * ch> C_qp;
    Eigen::Matrix<double, nc * ch, 1> lbC_qp;
    Eigen::Matrix<double, nc * ch, 1> ubC_qp;
    Eigen::Matrix<double, -1, -1> CI_, CE_;
    Eigen::VectorXd ci0_, ce0_;

    Vec12 F_;
    int tcounter;

    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;

    Eigen::Matrix<double, 3, 3> CrossProduct_A(Eigen::Matrix<double, 3, 1> A);
    Eigen::Matrix<double, 3, 3> Rz3(double theta);
    std::chrono::high_resolution_clock::time_point t1_prev;

    void setWeight();
    void solveQP();
    void Constraints();
};

#endif // MPC_H
