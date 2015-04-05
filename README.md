#SDL Tech Demo

install SDL with "sudo apt-get install libsdl2-2.0 libsdl2-dev"
then compile with "gcc -o test test.c `sdl2-config --cflags --libs`

The framerate looks like crap on the vm, but it's actually 60 frames per second.
We'll probably have to develop outside the vm for this reason.