SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: $(inter_tether_OBJS)

clean:
	rm -f $(inter_tether_OBJS) $(inter_tether_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_tether_OBJS:.o=.d)

-include $(inter_tether_OBJS:.o=.d)

