FROM rootproject/root
LABEL version=3.0.1

#-------------------------------------------------------------------------------
# Pre-configuration
#-------------------------------------------------------------------------------

# Necessary to get a nice log on the terminal
ENV TERM=xterm

#Install dependencies
RUN apt-get update && apt-get -y install git libomp-dev vim \
libgsl-dev rsync python3-pandas mayavi2 python3-matplotlib python3-scipy \
python3-tqdm zsh wget parallel

#Copy repository
RUN mkdir -p /home/MUSIC
COPY . /home/MUSIC

#Setup home folder
ENV HOME=/home
ENV USER=muisc
RUN mkdir -p /home/.parallel && touch /home/.parallel/will-cite

#Quality-of-life improvement: oh-my-zsh
WORKDIR /home
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.2/zsh-in-docker.sh)" -- \
    -t ys -p git

#-------------------------------------------------------------------------------
# Fix permissions
#-------------------------------------------------------------------------------
RUN chmod -R a+rw  /home
RUN chmod -R go-w /home/.oh-my-zsh

CMD zsh