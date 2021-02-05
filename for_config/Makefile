

include ../make.inc
code   = ../for_src
bin    = ../bin

all: 
	echo "do nothing"

%.o: %.f90
	cp $< ${code}/config.f90
	(cd ${code};make all )

%: %.o
	mv ${bin}/model.x ${bin}/${@}.x

clean: 
	(cd ${code}; make clean)


