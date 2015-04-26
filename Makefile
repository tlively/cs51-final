# This is the Makefile for building the project! yay!

all: game
	@echo "building all!"

game: game.c physics.o graphics.o
	gcc -Werror -o game game.c physics.o graphics.o `sdl2-config --cflags --libs` -lm

test: physics.o graphics.o
	gcc -Werror -o test test.c physics.o graphics.o `sdl2-config --cflags --libs` -lm && ./test

physics.o: physics.c
	gcc -Werror -o physics.o -c physics.c -lm

graphics.o: graphics.c
	gcc -Werror -o graphics.o -c graphics.c `sdl2-config --cflags --libs` -lm

clean:
	rm -f game test physics.o graphics.o

.PHONY: all clean test
