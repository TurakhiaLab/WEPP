import e, { useState as i } from "react";
import { r as O, i as q, j as R, o as b, d as $, P as S, e as M, f as P, g as w, h as p, k as G, l as v, I as z, m as H, n as C, M as J, G as E, u as y, B as g, p as x, L as V, b as W, c as K } from "./index-8ee3baa4.js";
import "react-dom";
var B = {}, Q = q;
Object.defineProperty(B, "__esModule", {
  value: !0
});
var I = B.default = void 0, U = Q(O()), X = R, Y = (0, U.default)(/* @__PURE__ */ (0, X.jsx)("path", {
  d: "M3 17.25V21h3.75L17.81 9.94l-3.75-3.75L3 17.25zM20.71 7.04c.39-.39.39-1.02 0-1.41l-2.34-2.34a.9959.9959 0 0 0-1.41 0l-1.83 1.83 3.75 3.75 1.83-1.83z"
}), "Create");
I = B.default = Y;
const Z = b(function({ rootModel: a, setIsAssemblyBeingEdited: c, setAssemblyBeingEdited: s }) {
  function t(n) {
    a.jbrowse.removeAssemblyConf(n);
  }
  const { assemblies: m } = a.jbrowse;
  return e.createElement(
    $,
    { component: S },
    e.createElement(
      M,
      null,
      e.createElement(
        P,
        null,
        e.createElement(
          w,
          null,
          e.createElement(p, null, "Name"),
          e.createElement(p, null, "Display name"),
          e.createElement(p, null, "Aliases"),
          e.createElement(p, null, "Actions")
        )
      ),
      e.createElement(G, null, m.map((n) => {
        const o = v.readConfObject(n, "name"), r = v.readConfObject(n, "displayName"), d = v.readConfObject(n, "aliases");
        return e.createElement(
          w,
          { key: o },
          e.createElement(p, null, o),
          e.createElement(p, null, r),
          e.createElement(p, null, d ? d.toString() : ""),
          e.createElement(
            p,
            null,
            e.createElement(
              z,
              { "data-testid": `${o}-edit`, onClick: () => {
                c(!0), s(n);
              } },
              e.createElement(I, { color: "primary" })
            ),
            e.createElement(
              z,
              { "data-testid": `${o}-delete`, onClick: () => t(o) },
              e.createElement(H, { color: "error" })
            )
          )
        );
      }))
    )
  );
}), ee = b(function({ adapterSelection: a, setAdapterSelection: c, adapterTypes: s }) {
  return e.createElement(C, { value: a, label: "Type", select: !0, helperText: "Type of adapter to use", fullWidth: !0, onChange: (t) => c(t.target.value) }, s.map((t) => e.createElement(J, { key: t, value: t }, t)));
}), te = b(({ adapterSelection: a, fastaLocation: c, setFastaLocation: s, faiLocation: t, setFaiLocation: m, gziLocation: n, setGziLocation: o, twoBitLocation: r, setTwoBitLocation: d, chromSizesLocation: f, setChromSizesLocation: L }) => a === "IndexedFastaAdapter" || a === "BgzipFastaAdapter" ? e.createElement(
  E,
  { container: !0, spacing: 2 },
  e.createElement(
    E,
    { item: !0 },
    e.createElement(y.FileSelector, { name: "fastaLocation", location: c, setLocation: (l) => s(l) })
  ),
  e.createElement(
    E,
    { item: !0 },
    e.createElement(y.FileSelector, { name: "faiLocation", location: t, setLocation: (l) => m(l) })
  ),
  a === "BgzipFastaAdapter" ? e.createElement(
    E,
    { item: !0 },
    e.createElement(y.FileSelector, { name: "gziLocation", location: n, setLocation: (l) => o(l) })
  ) : null
) : a === "TwoBitAdapter" ? e.createElement(
  E,
  { container: !0, spacing: 2 },
  e.createElement(
    E,
    { item: !0 },
    e.createElement(y.FileSelector, { name: "twoBitLocation", location: r, setLocation: (l) => d(l) })
  ),
  e.createElement(
    E,
    { item: !0 },
    e.createElement(y.FileSelector, { name: "chromSizesLocation (optional, can be added to speed up loading 2bit files with many contigs)", location: f, setLocation: (l) => L(l) })
  )
) : null), A = { uri: "" }, ae = b(function({ rootModel: a, setFormOpen: c }) {
  const s = [
    "IndexedFastaAdapter",
    "BgzipFastaAdapter",
    "TwoBitAdapter"
  ], [t, m] = i(""), [n, o] = i(""), [r, d] = i(s[0]), [f, L] = i(A), [l, k] = i(A), [h, j] = i(A), [F, _] = i(A), [T, N] = i(A);
  function D() {
    if (t === "")
      a.session.notify("Can't create an assembly without a name");
    else {
      c(!1);
      let u;
      r === "IndexedFastaAdapter" ? u = {
        name: t,
        displayName: n,
        sequence: {
          adapter: {
            type: "IndexedFastaAdapter",
            fastaLocation: f,
            faiLocation: l
          }
        }
      } : r === "BgzipFastaAdapter" ? u = {
        name: t,
        displayName: n,
        sequence: {
          adapter: {
            type: "BgzipFastaAdapter",
            fastaLocation: f,
            faiLocation: l,
            gziLocation: h
          }
        }
      } : r === "TwoBitAdapter" && (u = {
        name: t,
        displayName: n,
        sequence: {
          adapter: {
            type: "TwoBitAdapter",
            twoBitLocation: F,
            chromSizesLocation: T
          }
        }
      }), a.jbrowse.addAssemblyConf(u), a.session.notify(`Successfully added ${t} assembly to JBrowse 2`, "success");
    }
  }
  return e.createElement(
    "div",
    null,
    e.createElement(
      S,
      null,
      e.createElement(C, { id: "assembly-name", inputProps: { "data-testid": "assembly-name" }, label: "Assembly name", helperText: "The assembly name e.g. hg38", variant: "outlined", value: t, onChange: (u) => m(u.target.value) }),
      e.createElement(C, { id: "assembly-name", inputProps: { "data-testid": "assembly-display-name" }, label: "Assembly display name", helperText: 'A human readable display name for the assembly e.g. "Homo sapiens (hg38)"', variant: "outlined", value: n, onChange: (u) => o(u.target.value) }),
      e.createElement(ee, { adapterSelection: r, setAdapterSelection: d, adapterTypes: s }),
      e.createElement(te, { adapterSelection: r, fastaLocation: f, setFastaLocation: L, faiLocation: l, setFaiLocation: k, gziLocation: h, setGziLocation: j, twoBitLocation: F, setTwoBitLocation: _, chromSizesLocation: T, setChromSizesLocation: N })
    ),
    e.createElement(g, { variant: "contained", color: "secondary", startIcon: e.createElement(x, null), onClick: D }, "Create new assembly")
  );
}), ne = b(({ assembly: a }) => e.createElement(
  "div",
  { style: { maxHeight: 600, overflow: "auto" } },
  e.createElement(V, { model: { target: a } })
)), le = b(function({ rootModel: a, onClose: c }) {
  const [s, t] = i(!1), [m, n] = i(!1), [o, r] = i(), d = !s && !m;
  return e.createElement(
    y.Dialog,
    { open: !0, onClose: () => c(!1), title: "Assembly manager" },
    e.createElement(
      W,
      null,
      d ? e.createElement(Z, { rootModel: a, setIsAssemblyBeingEdited: n, setAssemblyBeingEdited: r }) : null,
      m ? e.createElement(ne, { assembly: o }) : null,
      s ? e.createElement(ae, { rootModel: a, setFormOpen: t }) : null
    ),
    e.createElement(
      K,
      null,
      s ? e.createElement(g, { variant: "contained", onClick: () => t(!1) }, "Back") : null,
      m ? e.createElement(g, { variant: "contained", onClick: () => n(!1) }, "Back") : null,
      d ? e.createElement(
        e.Fragment,
        null,
        e.createElement(g, { color: "secondary", variant: "contained", onClick: () => c(!1) }, "Close"),
        e.createElement(g, { variant: "contained", startIcon: e.createElement(x, null), onClick: () => t(!0) }, "Add new assembly")
      ) : null
    )
  );
}), ie = le;
export {
  ie as default
};
