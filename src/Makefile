all: LIME.so

LIME.so: LIME.cpp
	Rscript --vanilla -e "TMB::compile('LIME.cpp')"
    
clean:
	rm -rf *o