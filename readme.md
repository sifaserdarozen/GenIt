# GenIt
pseudo rtp payload generator experiment

### What is it?
Experiment code that uses Nvidia cuda in generating g711 / g722 rtp payloads. 
v0.1 version corresponds to code of 2015 Nvidia GPU technology conference poster
http://on-demand.gputechconf.com/gtc/2015/posters/GTC_2015_Signal___Audio_Processing_01_P5119_WEB.pdf

### How to build
In a Linux flavor, be sure to install *make*, *g++* and *nvidia-cuda-toolkit* package. Then clone or dowload code.
```
make
```

### How to run
Experiment may be performed through
```
./bin/genit
```

### License
Code of g722 encoder is highy resembles to ITU reference implementation that should be trated as IEFT wishes.
All remaining part is MIT.

