all: LIME.dll

LIME.dll: LIME.cpp
	Rscript --vanilla -e "TMB::compile('LIME.cpp')"

clean:
	rm -rf *o