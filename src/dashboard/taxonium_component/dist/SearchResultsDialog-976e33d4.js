import e from "react";
import { t as d, d as w, P as k, e as N, f as C, g as f, h as o, k as R, B as v, E as S, H as D, u as L, b as _, T as y, J as $, c as x } from "./index-8ee3baa4.js";
import "react-dom";
function I({ searchResults: u, assemblyName: E, model: a, handleClose: i }) {
  var r;
  const m = d.getSession(a), { pluginManager: b } = d.getEnv(m), { assemblyManager: h } = m, g = E || ((r = a.displayedRegions[0]) === null || r === void 0 ? void 0 : r.assemblyName), c = h.get(g);
  if (!c)
    throw new Error(`assembly ${g} not found`);
  if (!c.regions)
    throw new Error(`assembly ${g} regions not loaded`);
  function p(t) {
    var n;
    if (t) {
      const l = b.pluggableConfigSchemaType("track"), s = S(l, D(a), t);
      return ((n = s == null ? void 0 : s.name) === null || n === void 0 ? void 0 : n.value) || "";
    }
    return "";
  }
  async function T(t) {
    var n;
    try {
      const l = (n = c == null ? void 0 : c.regions) === null || n === void 0 ? void 0 : n.find((s) => t === s.refName);
      l ? (a.setDisplayedRegions([l]), a.showAllRegions()) : await a.navToLocString(t, g);
    } catch (l) {
      console.warn(l), m.notify(`${l}`, "warning");
    }
  }
  return e.createElement(
    w,
    { component: k },
    e.createElement(
      N,
      null,
      e.createElement(
        C,
        null,
        e.createElement(
          f,
          null,
          e.createElement(o, null, "Name"),
          e.createElement(o, { align: "right" }, "Location"),
          e.createElement(o, { align: "right" }, "Track"),
          e.createElement(o, { align: "right" })
        )
      ),
      e.createElement(R, null, u.map((t) => e.createElement(
        f,
        { key: `${t.getId()}` },
        e.createElement(o, { component: "th", scope: "row" }, t.getLabel()),
        e.createElement(o, { align: "right" }, t.getLocation()),
        e.createElement(o, { align: "right" }, p(t.getTrackId()) || "N/A"),
        e.createElement(
          o,
          { align: "right" },
          e.createElement(v, { onClick: async () => {
            try {
              const n = t.getLocation();
              if (n) {
                await T(n);
                const l = t.getTrackId();
                l && a.showTrack(l);
              }
            } catch (n) {
              console.error(n), m.notify(`${n}`, "error");
            }
            i();
          }, color: "primary", variant: "contained" }, "Go")
        )
      )))
    )
  );
}
function M({ model: u, assemblyName: E, searchQuery: a, searchResults: i, handleClose: r }) {
  return e.createElement(
    L.Dialog,
    { open: !0, maxWidth: "xl", onClose: r, title: "Search results" },
    e.createElement(_, null, i != null && i.length ? e.createElement(
      e.Fragment,
      null,
      e.createElement(
        y,
        null,
        "Showing results for ",
        e.createElement("b", null, a)
      ),
      e.createElement(I, { model: u, handleClose: r, assemblyName: E, searchResults: i })
    ) : e.createElement(
      y,
      null,
      "No results found for ",
      e.createElement("b", null, a)
    )),
    e.createElement($, null),
    e.createElement(
      x,
      null,
      e.createElement(v, { onClick: () => r(), color: "primary" }, "Cancel")
    )
  );
}
export {
  M as default
};
