# This is the Makefile for building the project! yay!

all: game
	@echo "building all!"

game: game.c physics.o graphics.o dynamic_array.o
	gcc -Werror -o game game.c physics.o graphics.o `sdl2-config --cflags --libs` -lm -std=c99

test: physics.o graphics.o
	gcc -Werror -o test test.c physics.o graphics.o `sdl2-config --cflags --libs` -lm -std=c99 && ./test

physics.o: physics.c dynamic_array.o
	gcc -Werror -o physics.o -c physics.c -lm -std=c99

graphics.o: graphics.c
	gcc -Werror -o graphics.o -c graphics.c `sdl2-config --cflags --libs` -lm -std=c99

dynamic_array.o: dynamic_array.c
	gcc -Werror -o dynamic_array.o -c dynamic_array.c `sdl2-config --cflags --libs` -lm -std=c99

clean:
	rm -f game test physics.o graphics.o

.PHONY: all clean test
