g='N=1e3, L=4, a=0, b=1, e=2, evt=1e-6, psd=1e-8, times=2e2'; d="run/m12"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst1($c, d=.30, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick

g='N=1e3, L=4, a=0, b=1, e=2, evt=1e-6, psd=1e-8, times=2e2'; d="run/m22"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst2($c, d=.20, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick

g='N=1e3, L=4, a=0, b=0, e=2, evt=1e-6, psd=1e-8, times=2e2'; d="run/m32"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst4($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst4($c, d=.20, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick

g='N=1e3, L=4, a=0, b=1, e=2, evt=1e-6, psd=1e-8, times=2e2'; d="run/m42"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst4($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst4($c, d=.20, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick

g='N=1e3, L=4, a=0, b=0, e=2, evt=1e-6, psd=1e-8, times=1e2'; d="run/m52"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst5($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst5($c, d=.10, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick

g='N=1e3, L=4, a=0, b=1, e=2, evt=1e-6, psd=1e-8, times=1e2'; d="run/m62"; rm -rf $d*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst5($c, d=.00, key=\"hd0\"); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst5($c, d=.10, key=\"hd2\"); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --ln 17q12.rds --tag ${d##*/} --par quick


for f in run/m??; do $f/sub.sh; done

for f in run/m??; do Rscript -e "source('rpt.R'); plt.pow('$f')"; done
