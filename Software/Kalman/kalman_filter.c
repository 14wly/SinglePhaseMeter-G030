/*
 * FileName : kalman_filter.c
 * Date     : 2022/3/31  
 */
 
#include "kalman_filter.h"
 
/*
 * @brief   
 *   初始化结构体kalman1_state
 *    我在这个init函数中设置了一些默认值:
 *    a = 1;
 *    h = 1;
 *
 *   NOTES: 请根据自己的要求更改 A,H,q,r 
 * @inputs  
 *   state - 卡尔曼结构体
 *   init_x -初始x状态值
 *   init_p -初始估计误差协方差
 * @outputs 
 * @retval  
 */
void kalman1_init(kalman1_state *state, float init_x, float init_p)
{
    state->x = init_x;
    state->p = init_p;
    state->A = 1;
    state->H = 1;
    state->q = 2e2;//10e-6;  /* 预测噪声协方差 */
    state->r = 5e2;//10e-5;  /* 测量误差协方差 */
}
 
/*
 * @brief   
 *   1维卡尔曼滤波
 * @inputs  
 *   state - 卡尔曼结构体
 *   z_measure - 测量值
 * @outputs 
 * @retval  
 *   估计结果
 */
float kalman1_filter(kalman1_state *state, float z_measure)
{
    /* 预测 */
    state->x = state->A * state->x;
    state->p = state->A * state->A * state->p + state->q;  /* p(n|n-1)=A^2*p(n-1|n-1)+q */
 
    /* 测量 */
    state->gain = state->p * state->H / (state->p * state->H * state->H + state->r);
    state->x = state->x + state->gain * (z_measure - state->H * state->x);
    state->p = (1 - state->gain * state->H) * state->p;
 
    return state->x;
}
 
/*
 * @brief   
 *   初始化结构体kalman1_state
 *   我在这个init函数中设置了一些默认值:
 *     A = {{1, 0.1}, {0, 1}};
 *     H = {1,0}; 
 *    
 *   NOTES: 请根据自己的要求更改 A,H,q,r 
 *
 * @inputs  
 * @outputs 
 * @retval  
 */
void kalman2_init(kalman2_state *state, float *init_x, float (*init_p)[2])
{
    state->x[0]    = init_x[0];
    state->x[1]    = init_x[1];
    state->p[0][0] = init_p[0][0];
    state->p[0][1] = init_p[0][1];
    state->p[1][0] = init_p[1][0];
    state->p[1][1] = init_p[1][1];
    //state->A       = {{1, 0.1}, {0, 1}};
    state->A[0][0] = 1;
    state->A[0][1] = 0.1;
    state->A[1][0] = 0;
    state->A[1][1] = 1;
    //state->H       = {1,0};
    state->H[0]    = 1;
    state->H[1]    = 0;
    //state->q       = {{10e-6,0}, {0,10e-6}};  /* 测量噪声协方差 */
    state->q[0]    = 10e-7;
    state->q[1]    = 10e-7;
    state->r       = 10e-7;  /* 估计误差协方差 */
}
 
/*
 * @brief   
 *   2维卡尔曼滤波
 * @inputs  
 *   state - 卡尔曼结构体
 *   z_measure - 测量值
 * @outputs 
 *   state->x[0] - 更新状态值，如角度、速度 
 *   state->x[1] - 更新状态值，如差角，加速度
 *   state->p - 更新估计误差约定矩阵
 * @retval  
 *   返回值等于state->x[0]，所以可能是角度或速度
 */
float kalman2_filter(kalman2_state *state, float z_measure)
{
    float temp0 = 0.0f;
    float temp1 = 0.0f;
    float temp = 0.0f;
 
    /* 预测 */
    state->x[0] = state->A[0][0] * state->x[0] + state->A[0][1] * state->x[1];
    state->x[1] = state->A[1][0] * state->x[0] + state->A[1][1] * state->x[1];
    /* p(n|n-1)=A^2*p(n-1|n-1)+q */
    state->p[0][0] = state->A[0][0] * state->p[0][0] + state->A[0][1] * state->p[1][0] + state->q[0];
    state->p[0][1] = state->A[0][0] * state->p[0][1] + state->A[1][1] * state->p[1][1];
    state->p[1][0] = state->A[1][0] * state->p[0][0] + state->A[0][1] * state->p[1][0];
    state->p[1][1] = state->A[1][0] * state->p[0][1] + state->A[1][1] * state->p[1][1] + state->q[1];
 
    /* 测量 */
    /* gain = p * H^T * [r + H * p * H^T]^(-1), H^T 意思是转置. */
    temp0 = state->p[0][0] * state->H[0] + state->p[0][1] * state->H[1];
    temp1 = state->p[1][0] * state->H[0] + state->p[1][1] * state->H[1];
    temp  = state->r + state->H[0] * temp0 + state->H[1] * temp1;
    state->gain[0] = temp0 / temp;
    state->gain[1] = temp1 / temp;
    /* x(n|n) = x(n|n-1) + gain(n) * [z_measure - H(n)*x(n|n-1)]*/
    temp = state->H[0] * state->x[0] + state->H[1] * state->x[1];
    state->x[0] = state->x[0] + state->gain[0] * (z_measure - temp); 
    state->x[1] = state->x[1] + state->gain[1] * (z_measure - temp);
 
    /* 更新 p: p(n|n) = [I - gain * H] * p(n|n-1) */
    state->p[0][0] = (1 - state->gain[0] * state->H[0]) * state->p[0][0];
    state->p[0][1] = (1 - state->gain[0] * state->H[1]) * state->p[0][1];
    state->p[1][0] = (1 - state->gain[1] * state->H[0]) * state->p[1][0];
    state->p[1][1] = (1 - state->gain[1] * state->H[1]) * state->p[1][1];
 
    return state->x[0];
}
