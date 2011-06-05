CXX = g++
OBJS = newick_parser/seqUtil.o newick_parser/Newickform.o

mas: myTree.h ${OBJS}
	make -C newick_parser
	${CXX} mas.cpp -o mas.out ${OBJS}

test: mas 35.tree 113.tree
	./mas.out 35.tree 113.tree agreement35and113.tree 2> /dev/null
	
