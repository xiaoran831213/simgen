g='N=1500, L=1, b=1, e=3, evt=1e-6, psd=1e-8, times=2e2'; d="run/m12"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {1..9..2}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.0, d=.5, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.2, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.2, d=.5, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none

g='L=1, M=3, b=1, e=3, evt=1e-6, psd=1e-8, times=2e2'; d="run/n12"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for n in "50+200*"{1..7}; do
        c="$g, N=$n, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.0, d=.8, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.3, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, a=.3, d=.8, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q28 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none

g='N=1500, L=1, b=1, e=2, evt=1e-6, psd=1e-8, times=200'; d="run/m22"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {1..9..2}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.0, d=.5, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.2, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.2, d=.5, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none

g='L=1, M=3, b=1, e=2, evt=1e-6, psd=1e-8, times=200'; d="run/n22"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for n in "50+200*"{1..7}; do
        c="$g, N=$n, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.0, d=.8, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.3, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, a=.3, d=.8, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q28 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none

g='N=1500, L=1, b=1, e=2, evt=1e-6, psd=1e-8, times=200'; d="run/m32"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {1..9..2}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.0, d=.5, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.2, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.2, d=.5, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none

g='L=1, M=3, b=1, e=2, evt=1e-6, psd=1e-8, times=200'; d="run/n32"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for n in "50+200*"{1..7}; do
        c="$g, N=$n, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.0, d=.0, key=\"H00\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.0, d=.8, key=\"H01\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.3, d=.0, key=\"H10\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst3($c, a=.3, d=.8, key=\"H11\"); $s'"
    done
done | hpcwp - -d$d -q28 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick --log none


for f in run/m??; do $f/sub.sh; done
for f in run/n??; do $f/sub.sh; done

for f in run/m??; do Rscript -e "source('rpt.R'); plt.pow('$f', x='M', xlab='# of Traits')"; done
for f in run/n??; do Rscript -e "source('rpt.R'); plt.pow('$f', x='N', xlab='Sample size')"; done
