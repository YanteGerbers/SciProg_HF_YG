$(DEST)/aux_boys_vec.o: $(DEST)/xkind.mod

$(DEST)/binom_coeff.o: $(DEST)/xkind.mod

$(DEST)/carmom_deriv.o: $(DEST)/xkind.mod

$(DEST)/carmom_hbra.o: $(DEST)/xkind.mod

$(DEST)/carmom_hrr_ket.o: $(DEST)/xkind.mod

$(DEST)/carmom_moment.o: $(DEST)/xkind.mod

$(DEST)/const_contr_gto.o: $(DEST)/xkind.mod

$(DEST)/const_contr_ints.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_carmom.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_delta.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_gaupot.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_nucpot.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_odist.o: $(DEST)/xkind.mod

$(DEST)/contr_cgto_value.o: $(DEST)/xkind.mod

$(DEST)/contr_csgto_carmom.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_carmom.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_delta.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_gaupot.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_nucpot.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_odist.o: $(DEST)/xkind.mod

$(DEST)/contr_sgto_value.o: $(DEST)/xkind.mod

$(DEST)/delta_geom.o: $(DEST)/xkind.mod

$(DEST)/delta_hket.o: $(DEST)/xkind.mod

$(DEST)/delta_moment.o: $(DEST)/xkind.mod

$(DEST)/gaupot_geom.o: $(DEST)/xkind.mod

$(DEST)/geom_part_one.o: $(DEST)/xkind.mod

$(DEST)/geom_part_zero.o: $(DEST)/xkind.mod

$(DEST)/geom_total.o: $(DEST)/xkind.mod

$(DEST)/get_address_list.o: $(DEST)/xkind.mod

$(DEST)/hgto_to_cgto.o: $(DEST)/xkind.mod

$(DEST)/hgto_to_lcgto.o: $(DEST)/xkind.mod

$(DEST)/hgto_to_sgto.o: $(DEST)/xkind.mod

$(DEST)/london_mom_hgto.o: $(DEST)/xkind.mod


$(DEST)/norm_contr_cgto.o: $(DEST)/xkind.mod

$(DEST)/norm_contr_sgto.o: $(DEST)/xkind.mod

$(DEST)/nucpot_geom.o: $(DEST)/xkind.mod

$(DEST)/nucpot_hbra.o: $(DEST)/xkind.mod

$(DEST)/nucpot_hket.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_carmom.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_delta.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_gaupot.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_nucpot.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_odist.o: $(DEST)/xkind.mod

$(DEST)/prim_hgto_value.o: $(DEST)/xkind.mod

$(DEST)/reorder_ints.o: $(DEST)/xkind.mod

$(DEST)/shell_scatter.o: $(DEST)/xkind.mod

$(DEST)/sort_cents.o: $(DEST)/xkind.mod

$(DEST)/trace_ints.o: $(DEST)/xkind.mod

$(DEST)/dump_info.o: $(DEST)/xkind.mod


$(DEST)/gen1int.o: $(DEST)/gen1int_carmom.mod $(DEST)/gen1int_gaupot.mod
$(DEST)/gen1int.o: $(DEST)/gen1int_geom.mod $(DEST)/gen1int_nucpot.mod
$(DEST)/gen1int.o: $(DEST)/gen1int_onehamil.mod $(DEST)/london_ao.mod
$(DEST)/gen1int.o: $(DEST)/xkind.mod

$(DEST)/gen1int_carmom.o: $(DEST)/gen1int_geom.mod $(DEST)/london_ao.mod
$(DEST)/gen1int_carmom.o: $(DEST)/xkind.mod

$(DEST)/gen1int_gaupot.o: $(DEST)/gen1int_geom.mod $(DEST)/london_ao.mod
$(DEST)/gen1int_gaupot.o: $(DEST)/xkind.mod

$(DEST)/gen1int_geom.o: $(DEST)/london_ao.mod $(DEST)/xkind.mod

$(DEST)/gen1int_nucpot.o: $(DEST)/gen1int_geom.mod $(DEST)/london_ao.mod
$(DEST)/gen1int_nucpot.o: $(DEST)/xkind.mod

$(DEST)/gen1int_onehamil.o: $(DEST)/gen1int_carmom.mod $(DEST)/gen1int_geom.mod
$(DEST)/gen1int_onehamil.o: $(DEST)/gen1int_nucpot.mod $(DEST)/london_ao.mod
$(DEST)/gen1int_onehamil.o: $(DEST)/xkind.mod

$(DEST)/london_ao.o: $(DEST)/xkind.mod


$(DEST)/xtimer.o: $(DEST)/xkind.mod


$(DEST)/module_interest_eri.o: $(DEST)/module_interest_hrr.mod
$(DEST)/module_interest_eri.o: $(DEST)/module_interest_osr.mod


$(DEST)/module_interest_one.o: $(DEST)/module_interest_osr.mod



$(DEST)/compute_integrals.o: $(DEST)/ao_basis.mod $(DEST)/gen1int.mod
$(DEST)/compute_integrals.o: $(DEST)/module_interest_eri.mod
$(DEST)/compute_integrals.o: $(DEST)/molecular_structure.mod




$(DEST)/main.o: $(DEST)/ao_basis.mod $(DEST)/compute_integrals.mod
$(DEST)/main.o: $(DEST)/diagonalization.mod $(DEST)/molecular_structure.mod
$(DEST)/main.o: $(DEST)/scf_convergence.mod $(DEST)/system_setup.mod

$(DEST)/system_setup.o: $(DEST)/ao_basis.mod $(DEST)/molecular_structure.mod

$(DEST)/ao_basis.mod: $(DEST)/ao_basis.o
$(DEST)/compute_integrals.mod: $(DEST)/compute_integrals.o
$(DEST)/diagonalization.mod: $(DEST)/diagonalization.o
$(DEST)/gen1int.mod: $(DEST)/gen1int.o
$(DEST)/gen1int_carmom.mod: $(DEST)/gen1int_carmom.o
$(DEST)/gen1int_gaupot.mod: $(DEST)/gen1int_gaupot.o
$(DEST)/gen1int_geom.mod: $(DEST)/gen1int_geom.o
$(DEST)/gen1int_nucpot.mod: $(DEST)/gen1int_nucpot.o
$(DEST)/gen1int_onehamil.mod: $(DEST)/gen1int_onehamil.o
$(DEST)/london_ao.mod: $(DEST)/london_ao.o
$(DEST)/module_interest_eri.mod: $(DEST)/module_interest_eri.o
$(DEST)/module_interest_hrr.mod: $(DEST)/module_interest_hrr.o
$(DEST)/module_interest_osr.mod: $(DEST)/module_interest_osr.o
$(DEST)/molecular_structure.mod: $(DEST)/molecular_structure.o
$(DEST)/scf_convergence.mod: $(DEST)/scf_convergence.o
$(DEST)/system_setup.mod: $(DEST)/system_setup.o
$(DEST)/xkind.mod: $(DEST)/xkind.o
