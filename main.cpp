/*
 * @Author: Fater7 
 * @Date: 2019-07-03 18:32:42 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 20:44:41
 */
/*
 * 基于简单遗传算法的正函数极值计算框架
 */

#include "GeneSFlow.h"

double Fx(double x)
{
    return x * sin(10*M_PI*x) + 2;
}

int main()
{
    GeneSFlow myGSF = GeneSFlow(70, 400, GSFRange(-1, 3), 0.000001, Fx);
    myGSF.SetProbForCross(0.7);
    myGSF.SetProbForMtp(0.007);
    double res = myGSF.start();
    cout << res << " " << Fx(res) << endl;
}