FROM vnmd/qsmxt_7.1.0:20240809
LABEL maintainer="Lee Reid <lee.reid1@uqconnect.edu.au>"
LABEL version=1.0
LABEL description="QSM QT, with an additional skull stripping step and real/imaginary dicom conversion"


RUN ls /dev
RUN apt-get update && apt-get install -y bash
SHELL ["/bin/bash", "-c"]

RUN apt install -y software-properties-common && add-apt-repository -y 'ppa:deadsnakes/ppa' && apt-get update -y && apt-get upgrade -y && apt-get install dcm2niix wget git -y 
COPY . /opt/
WORKDIR /opt

RUN ls /bin/
RUN /opt/install.sh


ENTRYPOINT [ "/opt/process-qsm.sh" ]

