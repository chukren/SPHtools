ifort -c SHTOOLS.f90
ifort -c FFTW3.f90
ifort -c PlmBar.f90
ifort -c SHRead.f90
ifort -c MakeGridDH.f90
#ifort -c sph.f90
#ifort -o readmarstopo readmarstopo.f90 PlmBar.o SHRead.o sph.o SHTOOLS.o MakeGridDH.o FFTW3.o  -L/home/youyir/tools/fftw3/lib -I/home/youyir/tools/fftw3/include -lfftw3 -m64 -lm -O3
ifort -o readmarstopo readmarstopo.f90 PlmBar.o SHRead.o SHTOOLS.o MakeGridDH.o FFTW3.o  -L/home/youyir/tools/fftw3/lib -I/home/youyir/tools/fftw3/include -lfftw3 -m64 -lm -O3

ifort -o readmarsmoho readmarsmoho.f90 PlmBar.o SHRead.o SHTOOLS.o MakeGridDH.o FFTW3.o  -L/home/youyir/tools/fftw3/lib -I/home/youyir/tools/fftw3/include -lfftw3 -m64 -lm -O3

#===
ifort -c SHPowerSpectra.f90
