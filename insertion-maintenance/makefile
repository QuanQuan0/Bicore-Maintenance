CC=g++
CPPFLAGS=-I
LDFLAGS=-g
BOOST_ROOT=PATH_TO_BOOST_ROOT
DEPS = bigraph.h gephi.h paper.h preprocessing.h random.h utility.h dyn_rebuild.h
OBJ = bigraph.o gephi.o main.o paper.o preprocessing.o random.o utility.o dyn_rebuild.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)$(BOOST_ROOT)

bicore: $(OBJ)
	$(CC) -std=c++1y -pthread $(LDFLAGS) -O3 -o $@ $^ $(CPPFLAGS)$(BOOST_ROOT)
