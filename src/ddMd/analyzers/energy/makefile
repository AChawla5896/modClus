BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_analyzers_energy_OBJS)

clean:
	rm -f $(ddMd_analyzers_energy_OBJS) $(ddMd_analyzers_energy_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_analyzers_energy_OBJS:.o=.d)

-include $(ddMd_analyzers_energy_OBJS:.o=.d)

