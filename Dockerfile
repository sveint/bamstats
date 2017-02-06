FROM centos:6

RUN yum install -y wget sudo zlib-devel

RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && \
    yum install -y devtoolset-2-gcc devtoolset-2-binutils

RUN curl -sSf https://static.rust-lang.org/rustup.sh > /home/rustup.sh
RUN chmod +x /home/rustup.sh
RUN sh /home/rustup.sh --disable-sudo -y


VOLUME ["/source"]
WORKDIR /source
