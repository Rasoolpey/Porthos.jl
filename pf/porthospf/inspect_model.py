"""Dump the PowerFactory model to JSON: network, machines, controllers, load flow.

Read-only on the user's study case: it activates the configured study case itself (no
working copy), reads attributes and runs a load flow (which changes no input data).
The attribute lists below are what the Porthos case needs; attributes that do not exist in
a PowerFactory version come out as null.
"""

from .session import name_of, safe_get

ELEMENT_ATTRS = {
    "ElmTerm": ["uknom", "outserv", "iUsage"],
    "ElmLne": ["dline", "nlnum", "outserv"],
    "ElmTr2": ["ntrcn", "nntap", "outserv"],
    "ElmLod": ["plini", "qlini", "scale0", "outserv"],
    "ElmShnt": ["qcapn", "qrean", "gparac", "outserv"],
    "ElmSym": ["pgini", "qgini", "usetp", "ip_ctrl", "ngnum", "outserv", "av_mode", "i_mot"],
    "ElmGenstat": ["pgini", "qgini", "usetp", "ip_ctrl", "ngnum", "outserv", "sgn"],
    "ElmXnet": ["bustp", "usetp", "phiini", "outserv"],
}
TYPE_ATTRS = {
    "ElmLne": ["uline", "rline", "xline", "bline", "gline", "sline"],
    "ElmTr2": ["strn", "utrn_h", "utrn_l", "uktr", "pcutr", "uk0tr", "curmg", "pfe", "dutap",
               "phitr", "nntap0", "ntpmn", "ntpmx"],
    "ElmLod": ["systp", "kpu", "kqu", "aP", "bP", "cP", "kpu0", "kpu1", "aQ", "bQ", "cQ",
               "kqu0", "kqu1", "loddy", "t1", "Tpf", "Tqf"],
    "ElmSym": ["sgn", "ugn", "cosn", "h", "dpu", "rstr", "xl", "xd", "xq", "xds", "xqs",
               "xdss", "xqss", "tds", "tqs", "tdss", "tqss", "tds0", "tqs0", "tdss0", "tqss0",
               "iturbo", "model_inp", "satur", "sg10", "sg12", "rotor", "xrl", "xrlq",
               "itimec", "iopt_dpu"],
}


def _terminals(obj):
    out = []
    for attr in ("bus1", "bus2", "bushv", "buslv"):
        cub = safe_get(obj, attr)
        term = safe_get(cub, "cterm") if cub is not None else None
        if term is not None:
            out.append(name_of(term))
    return out


def _dsl_param_names(dsl):
    raw = safe_get(safe_get(dsl, "typ_id"), "sParams")
    if raw is None:
        return []
    items = raw if isinstance(raw, (list, tuple)) else [raw]
    return [p.strip() for s in items for p in str(s).split(",") if p.strip()]


def _plain(v):
    """JSON-safe copy of an attribute value (objects become their names)."""
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    return name_of(v)


def _element(obj, cls):
    rec = {"name": name_of(obj), "class": cls, "terminals": _terminals(obj)}
    rec["attrs"] = {a: _plain(safe_get(obj, a)) for a in ELEMENT_ATTRS.get(cls, [])}
    typ = safe_get(obj, "typ_id")
    if typ is not None:
        rec["type"] = name_of(typ)
        rec["type_attrs"] = {a: _plain(safe_get(typ, a)) for a in TYPE_ATTRS.get(cls, [])}
    return rec


def dump(app):
    """The model of the active study case as a JSON-ready dict."""
    model = {"networks": [{"name": name_of(n), "frnom": safe_get(n, "frnom")}
                          for n in app.GetCalcRelevantObjects("*.ElmNet")]}
    for cls in ("ElmTerm", "ElmLne", "ElmTr2", "ElmLod", "ElmShnt", "ElmSym", "ElmGenstat",
                "ElmXnet"):
        model[cls] = [_element(o, cls) for o in app.GetCalcRelevantObjects("*." + cls)]
    model["ElmComp"] = [{"name": name_of(c), "frame": name_of(safe_get(c, "typ_id")),
                         "outserv": safe_get(c, "outserv"),
                         "slots": [name_of(s) for s in (safe_get(c, "pelm") or [])]}
                        for c in app.GetCalcRelevantObjects("*.ElmComp")]
    model["ElmDsl"] = [{"name": name_of(d), "type": name_of(safe_get(d, "typ_id")),
                        "outserv": safe_get(d, "outserv"),
                        "param_names": _dsl_param_names(d),
                        "params": _plain(safe_get(d, "params"))}
                       for d in app.GetCalcRelevantObjects("*.ElmDsl")]

    ldf = app.GetFromStudyCase("ComLdf")
    rc = ldf.Execute() if ldf is not None else None
    flow = {"rc": rc, "buses": [], "machines": []}
    if rc == 0:
        flow["buses"] = [{"name": name_of(t), "u": safe_get(t, "m:u"),
                          "phiu": safe_get(t, "m:phiu"), "U": safe_get(t, "m:U")}
                         for t in app.GetCalcRelevantObjects("*.ElmTerm")]
        for cls in ("ElmSym", "ElmGenstat"):
            flow["machines"] += [{"name": name_of(m), "class": cls,
                                  "P_MW": safe_get(m, "m:Psum:bus1"),
                                  "Q_Mvar": safe_get(m, "m:Qsum:bus1"),
                                  "u": safe_get(m, "m:u:bus1")}
                                 for m in app.GetCalcRelevantObjects("*." + cls)]
    model["load_flow"] = flow
    return model
