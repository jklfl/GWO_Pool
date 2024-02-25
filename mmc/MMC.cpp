#include <iostream>
#include <cmath>
#include <cstdlib>
const int arr_num = 1000000;
const int services = 5;
double expntl(double L) {
    double u = static_cast<double>(rand()) / RAND_MAX;
    return -L * log(u);
}
double MMC(double lam) {
    srand(time(NULL));
    const double Total_time = 20;
    const double lam = 4;
    const double mu = 3;
    const double arr_mean = 1 / lam;
    const double ser_mean = 1 / mu;
    double events[5][arr_num] = {0};
    double arr[arr_num];
    for (int i = 0; i < arr_num; ++i) {
        arr[i] = expntl(arr_mean);
    }
    for (int i = 0; i < arr_num; ++i) {
        events[0][i] = i > 0 ? events[0][i - 1] + arr[i] : arr[i];
        events[1][i] = expntl(ser_mean);
    }
    int numbers = 2;
    int member[] = {0, 1};
    events[2][0] = events[2][1] = 0;
    events[3][0] = events[0][0] + events[1][0];
    events[3][1] = events[0][1] + events[1][1];
    events[4][0] = events[4][1] = 1;
    double leave_time[] = {events[3][0], events[3][1]};
    int ser_num = 0;
    for (int k = 2; k < arr_num; ++k) {
        if (events[0][k] >= Total_time) {
            ser_num = k;
            break;
        } else {
            if (events[0][k] <= fmin(leave_time[0], leave_time[1])) {
                events[2][k] = fmin(leave_time[0], leave_time[1]) - events[0][k];
                events[3][k] = events[0][k] + events[1][k] + events[2][k];
            } else {
                events[2][k] = 0;
                events[3][k] = events[0][k] + events[1][k] + events[2][k];
            }
            leave_time[0] = fmax(leave_time[0], leave_time[1]);
            leave_time[1] = events[3][k];
            ++numbers;
            events[4][k] = 1;
        }
    }
    double c = lam / mu;
    double p = lam / (2 * mu);
    double p0 = 1 / (1 + c + c * c / (2 * (1 - p)));
    double Lq = pow(services * p, services) * p * p0 / (tgamma(services + 1) * pow(1 - p, 2));
    double Ls = Lq + c;
    double Ws = Ls / lam;
    double Wq = Lq / lam;
    return Lq;
}
