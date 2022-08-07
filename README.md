# FIR-Wiener-filtering-LMS-Algorithm
Exercise on Statistical Signal Processing course - Departement of Informatics and Telecommunications, National & Kapodistrian Uneversity of Athens.
# IMPLEMENTATION
- 1.mat file contains 300 samples of the signal d[n] and 300 samples of the signal x[n]=d[n]+u[n], where u[n] is white gaussian noise. FIR Wiener filtering was used to remove the noise and reconstruct the initial d[n] signal. The produced signal is depicted for different ranks of the filter (p=10,40 and 80), compared to d[n].

- 2.mat file contains 300 samples of the signal u[n] - white gaussian noise and y[n], which was produced after filtering u[n], using a filter with unknown parameters. To identify the unknown filter LMS algorithm was used. Different ranks of the filter were tested (p=2,3 and 4).

The report of the exercise, where the outputs are fully explained, is also uploaded on the repository in Greek.
