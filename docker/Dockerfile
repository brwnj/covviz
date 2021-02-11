FROM continuumio/miniconda3:4.6.14

RUN apt-get update && apt-get install -y procps
RUN conda create -n env -c conda-forge python=3.7 numpy>=1.16.2 jinja2>=2.10.1
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
RUN /opt/conda/envs/env/bin/pip install --no-cache-dir covviz==1.3.1
RUN wget -qO /usr/bin/goleft https://github.com/brentp/goleft/releases/download/v0.2.3/goleft_linux64
RUN chmod +x /usr/bin/goleft
