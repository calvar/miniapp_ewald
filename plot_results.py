import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)


def mean_stdev(a):
    av = a.sum()/len(a)
    stdev = np.sqrt( (a**2).sum()/len(a) - av**2)
    return av, stdev



#data--------------------------------
serial_real_time = np.array([99677.241, 99496.231, 100204.929])
serial_recip_time = np.array([62388.237, 62364.233, 62312.807])

omp_6_real_time = np.array([16784.935, 16786.759, 16780.025, 16623.18, 16755.747])
omp_6_recip_time = np.array([12642.276, 12670.332, 12631.381, 12620.815, 12613.547])

omp_14_real_time = np.array([7527.492, 7513.53, 7503.979, 7525.235, 7512.295])
omp_14_recip_time = np.array([5678.689, 5679.393, 5688.177, 5692.946, 5680.641])

omp_28_real_time = np.array([3830.935, 3814.398, 3840.55, 3830.026, 3825.146])
omp_28_recip_time = np.array([2958.292, 2971.852, 2970.49, 2959.804, 2979.072])

omp_64_real_time = np.array([2655.675, 2665.981, 2651.704, 2671.691, 2669.836])
omp_64_recip_time = np.array([2152.769, 2160.711, 2137.175, 2151.242, 2139.439])

cuda_64_real_time = np.array([314.152, 313.268, 309.155, 313.164, 312.806])
cuda_64_recip_time = np.array([149.696, 154.382, 154.42, 154.533, 133.201])

cuda_128_real_time = np.array([305.982, 309.276, 314.889, 314.305])
cuda_128_recip_time = np.array([142.543, 151, 152.997, 145.86])

cuda_256_real_time = np.array([306.877, 305.358, 313.742, 307.103, 308.929])
cuda_256_recip_time = np.array([181.248, 181.053, 181.135, 152.599, 159.609])

cuda_512_real_time = np.array([299.986, 310.68, 302.204, 300.884, 297.953])
cuda_512_recip_time = np.array([227.542, 231.943, 251.154, 242.501, 221.995])

#Compute average
ser_R_t, ser_R_std = mean_stdev(serial_real_time/1000)
ser_K_t, ser_K_std = mean_stdev(serial_recip_time/1000)

omp6_R_t, omp6_R_std = mean_stdev(omp_6_real_time/1000)
omp6_K_t, omp6_K_std = mean_stdev(omp_6_recip_time/1000)
omp14_R_t, omp14_R_std = mean_stdev(omp_14_real_time/1000)
omp14_K_t, omp14_K_std = mean_stdev(omp_14_recip_time/1000)
omp28_R_t, omp28_R_std = mean_stdev(omp_28_real_time/1000)
omp28_K_t, omp28_K_std = mean_stdev(omp_28_recip_time/1000)
omp64_R_t, omp64_R_std = mean_stdev(omp_64_real_time/1000)
omp64_K_t, omp64_K_std = mean_stdev(omp_64_recip_time/1000)

cuda64_R_t, cuda64_R_std = mean_stdev(cuda_64_real_time/1000)
cuda64_K_t, cuda64_K_std = mean_stdev(cuda_64_recip_time/1000)
cuda128_R_t, cuda128_R_std = mean_stdev(cuda_128_real_time/1000)
cuda128_K_t, cuda128_K_std = mean_stdev(cuda_128_recip_time/1000)
cuda256_R_t, cuda256_R_std = mean_stdev(cuda_256_real_time/1000)
cuda256_K_t, cuda256_K_std = mean_stdev(cuda_256_recip_time/1000)
cuda512_R_t, cuda512_R_std = mean_stdev(cuda_512_real_time/1000)
cuda512_K_t, cuda512_K_std = mean_stdev(cuda_512_recip_time/1000)

#Speedup
omp6_spdupR_t, omp6_spdupR_std = mean_stdev(ser_R_t/(omp_6_real_time/1000))
omp14_spdupR_t, omp14_spdupR_std = mean_stdev(ser_R_t/(omp_14_real_time/1000))
omp28_spdupR_t, omp28_spdupR_std = mean_stdev(ser_R_t/(omp_28_real_time/1000))
omp64_spdupR_t, omp64_spdupR_std = mean_stdev(ser_R_t/(omp_64_real_time/1000))
omp6_spdupK_t, omp6_spdupK_std = mean_stdev(ser_K_t/(omp_6_recip_time/1000))
omp14_spdupK_t, omp14_spdupK_std = mean_stdev(ser_K_t/(omp_14_recip_time/1000))
omp28_spdupK_t, omp28_spdupK_std = mean_stdev(ser_K_t/(omp_28_recip_time/1000))
omp64_spdupK_t, omp64_spdupK_std = mean_stdev(ser_K_t/(omp_64_recip_time/1000))


cuda64_spdupR_t, cuda64_spdupR_std = mean_stdev(omp64_R_t/(cuda_64_real_time/1000))
cuda128_spdupR_t, cuda128_spdupR_std = mean_stdev(omp64_R_t/(cuda_128_real_time/1000))
cuda256_spdupR_t, cuda256_spdupR_std = mean_stdev(omp64_R_t/(cuda_256_real_time/1000))
cuda512_spdupR_t, cuda512_spdupR_std = mean_stdev(omp64_R_t/(cuda_512_real_time/1000))
cuda64_spdupK_t, cuda64_spdupK_std = mean_stdev(omp64_K_t/(cuda_64_recip_time/1000))
cuda128_spdupK_t, cuda128_spdupK_std = mean_stdev(omp64_K_t/(cuda_128_recip_time/1000))
cuda256_spdupK_t, cuda256_spdupK_std = mean_stdev(omp64_K_t/(cuda_256_recip_time/1000))
cuda512_spdupK_t, cuda512_spdupK_std = mean_stdev(omp64_K_t/(cuda_512_recip_time/1000))

#print("{0:.3f}+/-{1:.3f}".format(cuda64_spdupR_t, cuda64_spdupR_std))

#plot arrays
P = [1, 6, 14, 28, 64]
Thr = [64, 128, 256, 512]
omp_time = [ser_R_t, omp6_R_t, omp14_R_t, omp28_R_t, omp64_R_t]
omp_time_sd = [ser_R_std, omp6_R_std, omp14_R_std, omp28_R_std, omp64_R_std]
omp_timek = [ser_K_t, omp6_K_t, omp14_K_t, omp28_K_t, omp64_K_t]
omp_timek_sd = [ser_K_std, omp6_K_std, omp14_K_std, omp28_K_std, omp64_K_std]
cuda_time = [cuda64_R_t, cuda128_R_t, cuda256_R_t, cuda512_R_t]
cuda_time_sd = [cuda64_R_std, cuda128_R_std, cuda256_R_std, cuda512_R_std]
cuda_timek = [cuda64_K_t, cuda128_K_t, cuda256_K_t, cuda512_K_t]
cuda_timek_sd = [cuda64_K_std, cuda128_K_std, cuda256_K_std, cuda512_K_std]

omp_spdup = [1, omp6_spdupR_t, omp14_spdupR_t, omp28_spdupR_t, omp64_spdupR_t]
omp_spdup_sd = [0, omp6_spdupR_std, omp14_spdupR_std, omp28_spdupR_std, omp64_spdupR_std]
omp_spdupk = [1, omp6_spdupK_t, omp14_spdupK_t, omp28_spdupK_t, omp64_spdupK_t]
omp_spdupk_sd = [0, omp6_spdupK_std, omp14_spdupK_std, omp28_spdupK_std, omp64_spdupK_std]
cuda_spdup = [cuda64_spdupR_t, cuda128_spdupR_t, cuda256_spdupR_t, cuda512_spdupR_t]
cuda_spdup_sd = [cuda64_spdupR_std, cuda128_spdupR_std, cuda256_spdupR_std, cuda512_spdupR_std]
cuda_spdupk = [cuda64_spdupK_t, cuda128_spdupK_t, cuda256_spdupK_t, cuda512_spdupK_t]
cuda_spdupk_sd = [cuda64_spdupK_std, cuda128_spdupK_std, cuda256_spdupK_std, cuda512_spdupK_std]

#Plot------------------------
fig1 = plt.figure(figsize=(25,20), dpi=30)
ax1 = fig1.add_subplot(111)

ax1.errorbar(P,omp_time,yerr=omp_time_sd,marker='o',label='Real')
ax1.errorbar(P,omp_timek,yerr=omp_timek_sd,marker='o',label='Reciprocal')
ax1.set_xlabel(r"$P$",fontsize=60)
ax1.set_ylabel(r"$\left<t\right>(sec.)$",fontsize=60)
ax1.tick_params(axis='both', which='major', labelsize=40)
ax1.legend(loc="best",fontsize=50,borderaxespad=0,bbox_to_anchor=(0.,0.,1.0,1.0))

plt.savefig("omp_time.png")
plt.draw()

##
fig2 = plt.figure(figsize=(25,20), dpi=30)
ax2 = fig2.add_subplot(111)

ax2.errorbar(Thr,cuda_time,yerr=cuda_time_sd,marker='o',label='Real')
ax2.errorbar(Thr,cuda_timek,yerr=cuda_timek_sd,marker='o',label='Reciprocal')
ax2.set_xlabel(r"threads/block",fontsize=60)
ax2.set_ylabel(r"$\left<t\right>(sec.)$",fontsize=60)
ax2.tick_params(axis='both', which='major', labelsize=40)
ax2.legend(loc="best",fontsize=50,borderaxespad=0,bbox_to_anchor=(0.,0.,1.0,1.0))

plt.savefig("cuda_time.png")
plt.draw()

##
fig3 = plt.figure(figsize=(25,20), dpi=30)
ax3 = fig3.add_subplot(111)

ax3.errorbar(P,omp_spdup,yerr=omp_spdup_sd,marker='o',label='Real')
ax3.errorbar(P,omp_spdupk,yerr=omp_spdupk_sd,marker='o',label='Reciprocal')
ax3.set_xlabel(r"$P$",fontsize=60)
ax3.set_ylabel(r"Speedup",fontsize=60)
ax3.tick_params(axis='both', which='major', labelsize=40)
ax3.legend(loc="best",fontsize=50,borderaxespad=0,bbox_to_anchor=(0.,0.,1.0,1.0))

plt.savefig("omp_speedup.png")
plt.draw()

##
fig4 = plt.figure(figsize=(25,20), dpi=30)
ax4 = fig4.add_subplot(111)

ax4.errorbar(Thr,cuda_spdup,yerr=cuda_spdup_sd,marker='o',label='Real')
ax4.errorbar(Thr,cuda_spdupk,yerr=cuda_spdupk_sd,marker='o',label='Reciprocal')
ax4.set_xlabel(r"Threads/block",fontsize=60)
ax4.set_ylabel(r"Speedup(wrt omp64)",fontsize=60)
ax4.tick_params(axis='both', which='major', labelsize=40)
ax4.legend(loc="best",fontsize=50,borderaxespad=0,bbox_to_anchor=(0.,0.,1.0,1.0))

plt.savefig("cuda_speedup.png")
plt.draw()
