# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)


__all__ = ["generate_roughness"]


def generate_roughness(roughness, roughness_ini, logger=logger):
    """ """
    roughness["frictionid"] = roughness.apply(
        lambda x: "%s_%s" % (x["frictiontype"], x["frictionvalue"]), axis=1
    )
    roughness_only = roughness[
        [roughness.index.name, "geometry", "frictiontype", "frictionvalue"]
    ]

    return roughness_only, roughness
