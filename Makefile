GRAPH = gnuplot

CC = clang
CFLAGS = -Wall -O0 -g 
LFLAGS = -O0 -g
LIBS = -lgsl -lm

star: WhiteSolver.o WhiteCode.o
	${CC} $(LFLAGS) -o $@ $^ $(LIBS)

data: star
	./star > data

M2Graph.png: WhiteGraph.gpl data
	$(GRAPH) WhiteGraph.gpl

clean : 
	rm -f *~
	rm -f *.o
	rm -f star

veryclean : clean
	rm -f data M2Graph.png
