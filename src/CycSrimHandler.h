#ifndef CYCSRIMHANDLER__H
#define CYCSRIMHANDLER__H

#include "CycSrim.h"

// To add new elements for consideration
// /opt/new/include/cycapp/CycSrim.h

// Reads the form "CycSrim::EMaterial" string from file and returns the corresponding
// material. Unfortunately, new elements have to be added here as needed. It's not
// efficient or clean, but it's what we're working with.
CycSrim::EMaterials evaluate_cycsrim_material(const std::string &material) {
    if(material == "CycSrim::SrimMaterialCarbon") {
        return CycSrim::SrimMaterialCarbon;
    } else if(material == "CycSrim::SrimMaterialAluminum") {
        return CycSrim::SrimMaterialAluminum;
    } else if(material == "CycSrim::SrimMaterialCD2") {
        return CycSrim::SrimMaterialCD2;
    } else if(material == "CycSrim::SrimMaterialSn") {
        return CycSrim::SrimMaterialSn;
    } else if(material == "CycSrim::SrimMaterialCu") {
        return CycSrim::SrimMaterialCu;
    } else if(material == "CycSrim::SrimMaterialAu") {
        return CycSrim::SrimMaterialAu;
    } else if(material == "CycSrim::SrimMaterialBoron") {
        return CycSrim::SrimMaterialBoron;
    } else {
        std::invalid_argument("CycSrimHandler::evaluate_cycsrim_material : invalid material " + material + ". Returning CycSrim::SrimMaterialAu");
        return CycSrim::SrimMaterialAu;
    }
}

/*
std::string return_cycsrim_material(int material) {
    if(material == 8) return "CycSrim::SrimMaterialCarbon";
    if(material == 4) return "CycSrim::SrimMaterialAluminum";
    if(material == 26) return "CycSrim::SrimMaterialCD2";
    if(material == 34) return "CycSrim::SrimMaterialSn";
    if(material == 50) return "CycSrim::SrimMaterialCu";
    if(material == 9) return "CycSrim::SrimMaterialAu";
    if(material == 51) return "CycSrim::SrimMaterialBoron";
    return "CycSrimHandler::return_cycsrim_material() : unknown material";
}
*/

#endif