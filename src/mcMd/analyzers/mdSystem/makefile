SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk

all: $(mcMd_analyzers_mdSystem_OBJS)

clean:
	rm -f $(mcMd_analyzers_mdSystem_OBJS) $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

-include $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

