/*
 * @Author: Fater7 
 * @Date: 2019-07-04 19:20:39 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 20:04:54
 */

#pragma once

typedef struct GSFRange
{
    GSFRange(double ts, double tl);
    GSFRange(): startValue(0), length(0){};
    double startValue;
    double length;
}GSFRange;