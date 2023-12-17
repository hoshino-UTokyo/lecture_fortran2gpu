subdirs=ACC_D_f_a  ACC_D_f_c  ACC_D_s_a  ACC_D_s_c ACC_U_f_a  ACC_U_f_c  ACC_U_s_a  ACC_U_s_c

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

