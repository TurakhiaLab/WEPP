import t, { useState as b, useReducer as ae, useMemo as F, useEffect as U } from "react";
import { G as W, n as le, a5 as re, I as T, bn as se, u as $, q as ie, ae as oe, T as ce, a9 as X, bo as de, p as me, aA as ue, o as Y, t as v, bp as fe, l as _, au as pe, bq as N, aR as he, X as ge, br as Ee, H as be, E as ve, aq as we, b as ke } from "./index-8ee3baa4.js";
import { u as ye, _ as Ce } from "./ResizeBar-ce571297.js";
import { S as Se, d as xe } from "./HierarchicalTrackSelector-a1670414.js";
import { d as J } from "./Clear-70a3606f.js";
import "react-dom";
function Te({ setFilterText: r, setUseShoppingCart: l, setHideSparse: u, setShowOptions: w, showOptions: g, hideSparse: f, useShoppingCart: p, filterText: E, model: s }) {
  const [i, m] = b(null);
  return t.createElement(
    t.Fragment,
    null,
    t.createElement(
      W,
      { container: !0, spacing: 4, alignItems: "center" },
      t.createElement(
        W,
        { item: !0 },
        t.createElement(le, { label: "Search...", value: E, onChange: (o) => r(o.target.value), InputProps: {
          endAdornment: t.createElement(
            re,
            { position: "end" },
            t.createElement(
              T,
              { onClick: () => r("") },
              t.createElement(J, null)
            )
          )
        } })
      ),
      t.createElement(
        W,
        { item: !0 },
        t.createElement(
          T,
          { onClick: (o) => m(o.currentTarget) },
          t.createElement(se, null)
        )
      ),
      t.createElement(
        W,
        { item: !0 },
        t.createElement(Se, { model: s })
      )
    ),
    t.createElement($.Menu, { anchorEl: i, open: !!i, onClose: () => m(null), onMenuItemClick: (o, h) => {
      h(), m(null);
    }, menuItems: [
      {
        label: "Add tracks to selection instead of turning them on/off",
        onClick: () => l(!p),
        type: "checkbox",
        checked: p
      },
      {
        label: "Hide sparse metadata columns",
        onClick: () => u(!f),
        checked: f,
        type: "checkbox"
      },
      {
        label: "Show extra table options",
        onClick: () => w(!g),
        checked: g,
        type: "checkbox"
      }
    ] })
  );
}
const Me = ie()((r) => ({
  facet: {
    margin: 0,
    marginLeft: r.spacing(2)
  },
  select: {
    marginBottom: r.spacing(2)
  }
}));
function Oe({ column: r, vals: l, width: u, dispatch: w, filters: g }) {
  const { classes: f } = Me(), [p, E] = b(!0);
  return t.createElement(
    oe,
    { key: r.field, className: f.facet, style: { width: u } },
    t.createElement(
      "div",
      { style: { display: "flex" } },
      t.createElement(ce, null, r.field),
      t.createElement(
        X,
        { title: "Clear selection on this facet filter" },
        t.createElement(
          T,
          { onClick: () => {
            w({ key: r.field, val: [] });
          }, size: "small" },
          t.createElement(J, null)
        )
      ),
      t.createElement(
        X,
        { title: "Minimize/expand this facet filter" },
        t.createElement(T, { onClick: () => E(!p), size: "small" }, p ? t.createElement(de, null) : t.createElement(me, null))
      )
    ),
    p ? t.createElement(ue, { multiple: !0, native: !0, className: f.select, value: g[r.field], onChange: (s) => {
      const { options: i } = s.target, m = [], o = i.length;
      for (let h = 0; h < o; h++)
        i[h].selected && m.push(i[h].value);
      w({ key: r.field, val: m });
    } }, l.sort((s, i) => s[0].localeCompare(i[0])).map(([s, i]) => t.createElement(
      "option",
      { key: s, value: s },
      s,
      " (",
      i,
      ")"
    ))) : null
  );
}
function Fe({ rows: r, columns: l, dispatch: u, filters: w, width: g }) {
  const f = l.slice(1), p = f.map(() => /* @__PURE__ */ new Map());
  for (const E of r)
    for (const [s, i] of f.entries()) {
      const m = p[s], o = `${E[i.field] || ""}`, h = m.get(o);
      o && (h === void 0 ? m.set(o, 1) : m.set(o, h + 1));
    }
  return t.createElement("div", null, f.map((E, s) => t.createElement(Oe, { key: E.field, vals: [...p[s]], column: E, width: g, dispatch: u, filters: w })));
}
function We(r) {
  return Object.entries(r).map(([l, u]) => typeof u == "string" ? l : "").filter((l) => !!l);
}
const _e = ["category", "adapter", "description"], j = 0.75, je = Y(function({ model: l }) {
  var u;
  const { assemblyNames: w, view: g, selection: f } = l, { pluginManager: p } = v.getEnv(l), { ref: E, scrollLeft: s } = ye(), [i, m] = b(""), [o, h] = b(!1), [k, H] = b(), [I, Q] = b(!1), [y, Z] = b(!0), [M, ee] = b(400), A = w[0], x = v.getSession(l), z = v.useDebounce(i, 400), [q, te] = ae((e, n) => ({ ...e, [n.key]: n.val }), {}), d = F(() => l.trackConfigurations(A).filter((e) => fe(z, e, x)).map((e) => {
    var n, a;
    const c = _.readConfObject(e, "metadata");
    return {
      id: e.trackId,
      conf: e,
      name: pe.getTrackName(e, x),
      category: (n = _.readConfObject(e, "category")) === null || n === void 0 ? void 0 : n.join(", "),
      adapter: (a = _.readConfObject(e, "adapter")) === null || a === void 0 ? void 0 : a.type,
      description: _.readConfObject(e, "description"),
      metadata: c,
      ...c
    };
  }), [A, l, z, x]), C = F(() => _e.filter((e) => y ? d.map((n) => n[e]).filter((n) => !!n).length > 5 : !0), [y, d]), S = F(() => [...new Set(d.flatMap((e) => We(e.metadata)))].filter((e) => y ? d.map((n) => n.metadata[e]).filter((n) => !!n).length > 5 : !0), [y, d]), D = F(() => ["name", ...C, ...S], [C, S]), [R, B] = b({
    name: v.measureGridWidth(d.map((e) => e.name), { maxWidth: 500, stripHTML: !0 }) + 15,
    ...Object.fromEntries(C.map((e) => [
      e,
      v.measureGridWidth(d.map((n) => n[e]), { maxWidth: 400, stripHTML: !0 })
    ])),
    ...Object.fromEntries(S.map((e) => [
      e,
      v.measureGridWidth(d.map((n) => n.metadata[e]), { maxWidth: 400, stripHTML: !0 })
    ]))
  }), [O, K] = b(Object.fromEntries(D.map((e) => [e, !0])));
  U(() => {
    K((e) => ({
      ...Object.fromEntries(D.map((n) => [n, !0])),
      ...e
    }));
  }, [D]), U(() => {
    B((e) => ({
      name: e.name,
      ...Object.fromEntries(C.filter((n) => O[n]).map((n) => [
        n,
        v.measureGridWidth(d.map((a) => a[n]), { stripHTML: !0, maxWidth: 400 })
      ])),
      ...Object.fromEntries(S.filter((n) => O[n]).map((n) => [
        n,
        v.measureGridWidth(d.map((a) => a.metadata[n]), { stripHTML: !0, maxWidth: 400 })
      ]))
    }));
  }, [S, O, C, y, d]);
  const G = v.useDebounce(R, 400), V = [
    {
      field: "name",
      hideable: !1,
      renderCell: (e) => {
        const { value: n, id: a, row: c } = e;
        return t.createElement(
          t.Fragment,
          null,
          t.createElement(N, { html: n }),
          t.createElement(
            T,
            { onClick: (L) => H({
              target: L.currentTarget,
              id: a,
              conf: c.conf
            }) },
            t.createElement(xe, null)
          )
        );
      },
      width: G.name || 100
      // can be undefined before useEffect update
    },
    ...C.map((e) => ({
      field: e,
      width: G[e] || 100,
      renderCell: (n) => {
        const { value: a } = n;
        return a ? t.createElement(N, { html: a }) : "";
      }
    })),
    ...S.map((e) => ({
      field: e,
      width: G[e] || 100,
      renderCell: (n) => {
        const { value: a } = n;
        return a ? t.createElement(N, { html: a }) : "";
      }
    }))
  ], P = g.tracks.map(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (e) => e.configuration.trackId
  ), ne = Object.entries(q).filter((e) => e[1].length > 0);
  return t.createElement(
    t.Fragment,
    null,
    k ? t.createElement(he, { anchorEl: k == null ? void 0 : k.target, menuItems: ((u = x.getTrackActionMenuItems) === null || u === void 0 ? void 0 : u.call(x, k.conf)) || [], onMenuItemClick: (e, n) => {
      n(), H(void 0);
    }, open: !!k, onClose: () => H(void 0) }) : null,
    t.createElement(Te, { setHideSparse: Z, setShowOptions: h, setFilterText: m, setUseShoppingCart: Q, hideSparse: y, showOptions: o, filterText: i, useShoppingCart: I, model: l }),
    t.createElement(
      "div",
      { ref: E, style: {
        display: "flex",
        overflow: "hidden",
        height: window.innerHeight * j,
        width: window.innerWidth * j
      } },
      t.createElement(
        "div",
        { style: {
          height: window.innerHeight * j,
          width: window.innerWidth * j - M
        } },
        t.createElement(Ce, { checkbox: !0, widths: Object.values(R), setWidths: (e) => B(Object.fromEntries(Object.entries(R).map((n, a) => [
          n[0],
          e[a]
        ]))), scrollLeft: s }),
        t.createElement(ge, { rows: d.filter((e) => ne.every(([n, a]) => a.includes(e[n]))), columnVisibilityModel: O, onColumnVisibilityModelChange: (e) => K(e), columnHeaderHeight: 35, checkboxSelection: !0, disableRowSelectionOnClick: !0, keepNonExistentRowsSelected: !0, onRowSelectionModelChange: (e) => {
          if (I) {
            const n = be(l), a = p.pluggableConfigSchemaType("track"), c = e.map((L) => ve(a, n, L));
            l.setSelection(c);
          } else {
            const n = P, a = e;
            Ee(() => {
              n.filter((c) => !a.includes(c)).map((c) => g.hideTrack(c)), a.filter((c) => !n.includes(c)).map((c) => g.showTrack(c));
            });
          }
        }, rowSelectionModel: I ? f.map((e) => e.trackId) : P, slots: { toolbar: o ? we : null }, slotProps: {
          toolbar: { printOptions: { disableToolbarButton: !0 } }
        }, columns: V, rowHeight: 25 })
      ),
      t.createElement($.ResizeHandle, { vertical: !0, onDrag: (e) => ee(M - e), style: { background: "grey", width: 5 } }),
      t.createElement(
        "div",
        { style: { width: M, overflowY: "auto", overflowX: "hidden" } },
        t.createElement(Fe, { width: M - 10, rows: d, columns: V, dispatch: te, filters: q })
      )
    )
  );
});
function He(r) {
  const { handleClose: l } = r;
  return t.createElement(
    $.Dialog,
    { open: !0, onClose: l, maxWidth: "xl", title: "Faceted track selector" },
    t.createElement(
      ke,
      null,
      t.createElement(je, { ...r })
    )
  );
}
const $e = Y(He);
export {
  $e as default
};
