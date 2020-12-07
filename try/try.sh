g='N=1e3, a=0, b=.1, e=1, psd=1e-4, times=2e2'; d="run/a04"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1100}; do
    for m in {1..5}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=0, evt=1e-4, key=\"YT0\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=0, evt=1e-8, key=\"YT0\", tag=\"PD8\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=1, evt=1e-4, key=\"YT1\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=1, evt=1e-8, key=\"YT1\", tag=\"PD8\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}

g='N=1e3, M=4, a=0, e=1, psd=1e-4, times=2e2'; d="run/b00"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1100}; do
    for b in 0.{1..5}; do
        c="$g, b=$b, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=0, evt=1e-4, key=\"YT0\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=0, evt=1e-8, key=\"YT0\", tag=\"PD8\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=1, evt=1e-4, key=\"YT1\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst1($c, ydt=1, evt=1e-8, key=\"YT1\", tag=\"PD8\"); $s'"
    done
done | hpcwp - -d$d -q20 -m8 -p4 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}

# ---------------- collider test ---------------- #
g='N=1e3, es1=1, es2=1, psd=1e-4, times=2e2'; d="run/c04"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1100}; do
    for m in {1..5}; do
        c="$g, ydm=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=0, evt=1e-4, key=\"YT0\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=0, evt=1e-8, key=\"YT0\", tag=\"PD8\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=1, evt=1e-4, key=\"YT1\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=1, evt=1e-8, key=\"YT1\", tag=\"PD8\"); $s'"
    done
done | hpcwp - -d$d -q30 -m8 -p4 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}

g='N=1e3, es1=1, es2=1, psd=1e-8, times=2e2'; d="run/c08"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1000..1100}; do
    for m in {1..5}; do
        c="$g, ydm=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=0, evt=1e-4, key=\"YT0\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=0, evt=1e-8, key=\"YT0\", tag=\"PD8\"); $s'"

	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=1, evt=1e-4, key=\"YT1\", tag=\"PD4\"); $s'"
	echo "Rscript -e 'source(\"ini.R\"); r=tst2($c, ydt=1, evt=1e-8, key=\"YT1\", tag=\"PD8\"); $s'"
    done
done | hpcwp - -d$d -q30 -m8 -p4 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}

# ---------------- mvQTL with interation ---------------- #
g='N=1e3, a=0, b=1, e=2, times=2e2'; d="run/mx6"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst6($c, d=0, evt=1e-6, psd=1e-7); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst6($c, d=1, evt=1e-6, psd=1e-7); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p4 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}

g='N=1e3, b=1, d=0, e=2, times=2e2'; d="run/mx7"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1001..1100}; do
    for m in {2..6}; do
        c="$g, M=$m, seed=$i"; s="saveRDS(r, \"{n:05d}.rds\")"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst6($c, a=.00, evt=1e-5, psd=1e-6); $s'"
	echo "time Rscript -e 'source(\"ini.R\"); r=tst6($c, a=.05, evt=1e-5, psd=1e-6); $s'"
    done
done | hpcwp - -d$d -q10 -m8 -p1 --wtm 1 --cp "*.R" --cp R --tag ${d##*/}
