mcMd_user_=\
    mcMd/user/Environment.cpp\
    mcMd/user/hotClusterIdentifier.cpp\
    mcMd/user/hotBeadClustering.cpp

mcMd_user_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_user_))
mcMd_user_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_user_:.cpp=.o))

