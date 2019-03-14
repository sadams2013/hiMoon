FROM python:3.7

RUN apt-get update && apt-get -y upgrade && \
    apt install -y python3-pip samtools && \
    pip install pysam numba 

COPY setup.py config.ini scripts /hiMoonPGx/

COPY hiMoon hiMoonPGx/hiMoon

RUN cd hiMoonPGx && python setup.py install

CMD ["bash"]