import e, { useState as g } from "react";
import { t as C, Y as y, T as w, s as F, C as S, X as A, aq as H, n as q, q as I, o as L, u as R, b as x, c as M, B, v as j, ar as G, Z as _, P as J, J as N, as as Q } from "./index-8ee3baa4.js";
import { u as W, _ as z } from "./ResizeBar-ce571297.js";
import "react-dom";
function K({ columns: o, filter: n, setFilter: l }) {
  return e.createElement(
    e.Fragment,
    null,
    e.createElement(w, null, "These filters can use a plain text search or regex style query, e.g. in the genotype field, entering 1 will query for all genotypes that include the first alternate allele e.g. 0|1 or 1|1, entering [1-9]\\d* will find any non-zero allele e.g. 0|2 or 2/33"),
    o.map(({ field: a }) => e.createElement(q, { key: `filter-${a}`, placeholder: `Filter ${a}`, value: n[a] || "", onChange: (r) => l({ ...n, [a]: r.target.value }) }))
  );
}
function U(o) {
  var n;
  const { feature: l, descriptions: a = {} } = o, { ref: r, scrollLeft: i } = W(), [c, u] = g({}), d = l.samples || {}, m = Object.entries(d);
  let t, h = [];
  const s = Object.keys(c);
  try {
    h = m.map((p) => ({
      ...Object.fromEntries(Object.entries(p[1]).map((E) => [E[0], `${E[1]}`])),
      sample: p[0],
      id: p[0]
    })).filter((p) => s.length ? s.every((E) => {
      const k = c[E];
      return k ? p[E].match(new RegExp(k, "i")) : !0;
    }) : !0);
  } catch (p) {
    t = p;
  }
  const f = ["sample", ...Object.keys(((n = m[0]) === null || n === void 0 ? void 0 : n[1]) || {})], [v, O] = g(!1), [b, V] = g(f.map((p) => C.measureGridWidth(h.map((E) => E[p])))), T = f.map((p, E) => {
    var k, D;
    return {
      field: p,
      description: (D = (k = a.FORMAT) === null || k === void 0 ? void 0 : k[p]) === null || D === void 0 ? void 0 : D.Description,
      width: b[E]
    };
  });
  return m.length ? e.createElement(
    y.BaseCard,
    { ...o, title: "Samples" },
    t ? e.createElement(w, { color: "error" }, `${t}`) : null,
    e.createElement(F, { control: e.createElement(S, { checked: v, onChange: (p) => O(p.target.checked) }), label: e.createElement(w, { variant: "body2" }, "Show options") }),
    v ? e.createElement(K, { setFilter: u, columns: T, filter: c }) : null,
    e.createElement(
      "div",
      { ref: r },
      e.createElement(z, { widths: b, setWidths: V, scrollLeft: i }),
      e.createElement(A, { rows: h, hideFooter: h.length < 100, columns: T, disableRowSelectionOnClick: !0, rowHeight: 25, columnHeaderHeight: 35, disableColumnMenu: !0, slots: { toolbar: v ? H : null }, slotProps: {
        toolbar: { printOptions: { disableToolbarButton: !0 } }
      } })
    )
  ) : null;
}
function P({ rows: o, columns: n }) {
  const [i, c] = g(!1), u = Math.min(o.length, 100) * 25 + 100 + (i ? 50 : 0) + 50;
  return o.length ? e.createElement(
    "div",
    null,
    e.createElement(F, { control: e.createElement(S, { checked: i, onChange: (d) => c(d.target.checked) }), label: e.createElement(w, { variant: "body2" }, "Show options") }),
    e.createElement(
      "div",
      { style: {
        height: u,
        width: "100%"
      } },
      e.createElement(A, { rowHeight: 25, rows: o, columns: n, slots: { toolbar: i ? H : null } })
    )
  ) : null;
}
function X({ feature: o, descriptions: n }) {
  var l, a, r, i;
  const c = ((i = (r = (a = (l = n == null ? void 0 : n.INFO) === null || l === void 0 ? void 0 : l.CSQ) === null || a === void 0 ? void 0 : a.Description) === null || r === void 0 ? void 0 : r.match(/.*Format: (.*)/)) === null || i === void 0 ? void 0 : i[1].split("|")) || [], u = o.INFO.CSQ || [], d = u.map((t, h) => ({
    id: h,
    ...Object.fromEntries(t.split("|").map((s, f) => [c[f], s]))
  })) || [], m = c.map((t) => ({
    field: t
  }));
  return u.length ? e.createElement(
    y.BaseCard,
    { title: "CSQ table" },
    e.createElement(P, { rows: d, columns: m })
  ) : null;
}
function Y({ feature: o, descriptions: n }) {
  var l, a, r, i;
  const c = ((i = (r = (a = (l = n == null ? void 0 : n.INFO) === null || l === void 0 ? void 0 : l.ANN) === null || a === void 0 ? void 0 : a.Description) === null || r === void 0 ? void 0 : r.match(/.*Functional annotations:'(.*)'$/)) === null || i === void 0 ? void 0 : i[1].split("|")) || [], u = o.INFO.ANN || [], d = u.map((t, h) => ({
    id: h,
    ...Object.fromEntries(t.split("|").map((s, f) => [c[f], s]))
  })) || [], m = c.map((t) => ({
    field: t
  }));
  return u.length ? e.createElement(
    y.BaseCard,
    { title: "ANN table" },
    e.createElement(P, { rows: d, columns: m })
  ) : null;
}
const Z = I()({
  block: {
    display: "block"
  }
}), ee = L(function({ model: n, handleClose: l, feature: a, viewType: r }) {
  const { classes: i } = Z(), [c, u] = g(!0), [d, m] = g(!0);
  return e.createElement(
    R.Dialog,
    { open: !0, onClose: l, title: "Breakpoint split view options" },
    e.createElement(
      x,
      null,
      e.createElement(F, { className: i.block, control: e.createElement(S, { checked: c, onChange: (t) => u(t.target.checked) }), label: "Copy tracks into the new view" }),
      e.createElement(F, { className: i.block, control: e.createElement(S, { checked: d, onChange: (t) => m(t.target.checked) }), label: "Mirror tracks vertically in vertically stacked view" })
    ),
    e.createElement(
      M,
      null,
      e.createElement(B, { onClick: () => {
        const { view: t } = n, h = C.getSession(n);
        try {
          let f = function(O) {
            return O.map((b) => ({
              ...b,
              id: `${b.trackId}-${Math.random()}`
            }));
          };
          const s = r.snapshotFromBreakendFeature(a, t);
          s.views[0].offsetPx -= t.width / 2 + 100, s.views[1].offsetPx -= t.width / 2 + 100, s.featureData = a;
          const v = j(t.tracks);
          s.views[0].tracks = f(v), s.views[1].tracks = f(d ? [...v].reverse() : v), h.addView("BreakpointSplitView", s);
        } catch (s) {
          console.error(s), h.notify(`${s}`);
        }
        l();
      }, variant: "contained", color: "primary", autoFocus: !0 }, "OK"),
      e.createElement(B, { onClick: () => l(), color: "secondary", variant: "contained" }, "Cancel")
    )
  );
});
function $(o) {
  const { model: n, locStrings: l, feature: a } = o, r = C.getSession(n), { pluginManager: i } = C.getEnv(r), [c, u] = g(!1);
  let d;
  try {
    d = i.getViewType("BreakpointSplitView");
  } catch {
  }
  const m = new G(a);
  return e.createElement(
    y.BaseCard,
    { ...o, title: "Breakends" },
    e.createElement(w, null, "Link to linear view of breakend endpoints"),
    e.createElement("ul", null, l.map((t) => e.createElement(
      "li",
      { key: `${JSON.stringify(t)}` },
      e.createElement(_, { href: "#", onClick: (h) => {
        var s;
        h.preventDefault();
        const { view: f } = n;
        try {
          if (f)
            (s = f.navToLocString) === null || s === void 0 || s.call(f, t);
          else
            throw new Error("No view associated with this feature detail panel anymore");
        } catch (v) {
          console.error(v), r.notify(`${v}`);
        }
      } }, `LGV - ${t}`)
    ))),
    d ? e.createElement(
      "div",
      null,
      e.createElement(w, null, "Launch split views with breakend source and target"),
      e.createElement("ul", null, l.map((t) => e.createElement(
        "li",
        { key: `${JSON.stringify(t)}` },
        e.createElement(_, { href: "#", onClick: (h) => {
          h.preventDefault(), u(!0);
        } }, `${a.refName}:${a.start} // ${t} (split view)`)
      ))),
      c ? e.createElement(ee, { model: n, feature: m, viewType: d, handleClose: () => {
        u(!1);
      } }) : null
    ) : null
  );
}
function te(o) {
  const { model: n } = o, { featureData: l, descriptions: a } = n, r = JSON.parse(JSON.stringify(l)), { samples: i, ...c } = r, u = {
    CHROM: "chromosome: An identifier from the reference genome",
    POS: "position: The reference position, with the 1st base having position 1",
    ID: "identifier: Semi-colon separated list of unique identifiers where available",
    REF: "reference base(s): Each base must be one of A,C,G,T,N (case insensitive).",
    ALT: "alternate base(s): Comma-separated list of alternate non-reference alleles",
    QUAL: "quality: Phred-scaled quality score for the assertion made in ALT",
    FILTER: "filter status: PASS if this position has passed all filters, otherwise a semicolon-separated list of codes for filters that fail"
  };
  return e.createElement(
    J,
    { "data-testid": "variant-side-drawer" },
    e.createElement(y.FeatureDetails, { feature: c, descriptions: { ...u, ...a }, ...o }),
    e.createElement(N, null),
    e.createElement(X, { feature: c, descriptions: a }),
    e.createElement(N, null),
    e.createElement(Y, { feature: c, descriptions: a }),
    e.createElement(N, null),
    r.type === "breakend" ? e.createElement($, { feature: r, locStrings: r.ALT.map((d) => {
      var m;
      return ((m = Q(d)) === null || m === void 0 ? void 0 : m.MatePosition) || "";
    }), model: n }) : null,
    r.type === "translocation" ? e.createElement($, { feature: r, model: n, locStrings: [`${r.INFO.CHR2[0]}:${r.INFO.END}`] }) : null,
    e.createElement(U, { feature: r, ...o, descriptions: a })
  );
}
const oe = L(te);
export {
  oe as default
};
