#include "config.h"
#include "catch.hpp"

#include "guess.h"
#include "growth.h"
#include "management.h"

static const double TOLERANCE = 1e-10;
static const double COARSE_ROOT_OUTTAKE = 0.1;
static const double TWIG_OUTTAKE = 0.4;
static const double HARVEST_EFFICIENCY = 0.9;

static Individual create_individual(Pft &tree, Patch &patch, double dens, double biomass, double height) {
    Individual& indiv=patch.vegetation.createobj(tree, patch.vegetation);
    indiv.densindiv=dens;
    indiv.alive=true;
    double bminit=biomass;
    double ltor=4.0;
    allocation_init(bminit,ltor,indiv);
    allometry(indiv);
    indiv.height=height;
    return indiv;
}


TEST_CASE("Management harvests one tree", "[harvest][emission][litter]") {
    ifslowharvestpool = false;

    pftlist.killall();

    npft = 1;
    Pft tree = pftlist.createobj();
    tree.k_allom2 = 40.0;
    tree.k_allom3 = 0.67;
    tree.stem_frac = 0.65;
    tree.twig_frac = 0.13;
    tree.lifeform = TREE;

    Soiltype soilType = Soiltype();
    Gridcell gc = Gridcell();
    Stand stand = Stand(0, &gc, soilType, FOREST, 1);
    Patch &patch = stand.getobj();

    double initial_wood_biomass = 12.0;
    Individual indiv = create_individual(tree, patch, 1.0, initial_wood_biomass, 30.0);

    SECTION( "Harvest split up into litter and atmosphere" ) {

        harvest_wood(indiv, 1.0, HARVEST_EFFICIENCY, TWIG_OUTTAKE, COARSE_ROOT_OUTTAKE, false);


        // entire tree gone
        REQUIRE(indiv.cmass_wood() == Approx(0.0).margin(TOLERANCE));


        double expected_wood_harvest = initial_wood_biomass * 0.65 * HARVEST_EFFICIENCY;
        double residue_harvest =
                TWIG_OUTTAKE * 0.13 * initial_wood_biomass + COARSE_ROOT_OUTTAKE * 0.22 * initial_wood_biomass;

        REQUIRE(patch.fluxes.get_annual_flux(Fluxes::HARVESTC) ==
                Approx(expected_wood_harvest + residue_harvest).margin(TOLERANCE));

        // everything not harvested goes into litter
        REQUIRE(indiv.patchpft().total_litter() ==
                Approx(initial_wood_biomass - expected_wood_harvest - residue_harvest).margin(TOLERANCE));

        // no slow pool in this case
        REQUIRE(patch.pft[0].cmass_harvested_products_slow == Approx(0).margin(TOLERANCE));
    }


    SECTION( "Parts of the harvest go into the product pool" , "[slowpool]") {
        ifslowharvestpool = true;
        tree.harvest_slow_frac = 0.6;

        harvest_wood(indiv, 1.0, HARVEST_EFFICIENCY, TWIG_OUTTAKE, COARSE_ROOT_OUTTAKE, false);

        // entire tree gone
        REQUIRE(indiv.cmass_wood() == Approx(0.0).margin(TOLERANCE));

        double expected_wood_harvest = initial_wood_biomass * 0.65 * HARVEST_EFFICIENCY;
        double expected_slow_pool_harvest = expected_wood_harvest * tree.harvest_slow_frac;
        REQUIRE(indiv.patchpft().cmass_harvested_products_slow == Approx(expected_slow_pool_harvest).margin(TOLERANCE));

        double residue_harvest = TWIG_OUTTAKE * 0.13 * initial_wood_biomass + COARSE_ROOT_OUTTAKE * 0.22 * initial_wood_biomass;
        double expected_litter = initial_wood_biomass - expected_wood_harvest - residue_harvest;
        REQUIRE(indiv.patchpft().total_litter() == Approx(expected_litter).margin(TOLERANCE));

        // rest goes to litter
        double emissions = initial_wood_biomass - expected_litter - expected_slow_pool_harvest;
        REQUIRE(patch.fluxes.get_annual_flux(Fluxes::HARVESTC) == Approx(emissions).margin(TOLERANCE));
    }
}

