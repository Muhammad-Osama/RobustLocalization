# RobustLocalization
A method to localize passsive and auxiliary nodes in a wireless network in presence of outlier measurments. 

# Description of the .m files
There are three main scripts. Each of these scripts demonstrate the compares the proposed robust localization with the standard method in presence of non-line-of-sight measurements for a test scenario used for indoor localization called 3GPP:
1) ThreeGPP_tdoa_github.m : This script compares how the standard and robust methods when 'time-difference-of-arrival (tdoa)' measurements are made at the passive node. The performance is compared by estimating position for 50 different data draws over 50 Monte carlo simulations and the plotting the CDF of the absolute error in [m].
2) ThreeGPP_toa_github.m : This script compares how the standard and robust methods when 'time-of-arrival (tdoa)' measurements are made at the receiver node (which is also a transmitter in this case). The performance is compared by estimating position for 50 different data draws over 50 Monte carlo simulations and the plotting the CDF of the absolute error in [m].
3) ThreeGPP_sbp_github.m : This script compares how the standard and robust methods when anchor nodes transmitt in a schedule based localization scenario. The measurements are difference between time of arrival of different transmitters' signals at the passive node. This technique also allows the localization of auxiliary nodes. However, the script only demonstrates the localization of passive node. For details on schedule based localization, refer to 'Schedule-Based Localization in Asynchronous Wireless Networks' by Zachariah et al.

The rest of the scripts are supporting functions. Simply run the scripts to get results. 

For details regarding the robust localization method, please refer to our paper on arxiv: 'Robust Localization in Wireless Networks from Corrupted Signals', Muhammad Osama, Dave Zachariah and Peter Stoica.
