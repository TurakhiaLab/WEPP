import e, { useState as d, useEffect as T } from "react";
import { r as I, i as M, j as G, o as z, q as B, t as b, u as P, b as H, T as A, w as W, x as V, n as J, F as K, s as $, C as N, c as Q, B as x, y as U, z as X, A as Y, l as Z } from "./index-8ee3baa4.js";
import "react-dom";
var g = {};
Object.defineProperty(g, "__esModule", { value: !0 });
var j = g.formatSeqFasta = g.formatFastaLines = void 0;
function D(n) {
  return n.replace(/(.{1,80})/g, `$1
`).trimEnd();
}
g.formatFastaLines = D;
function ee(n) {
  return n.map((s) => `>${s.header}
${D(s.seq)}`).join(`
`);
}
j = g.formatSeqFasta = ee;
var F = {}, te = M;
Object.defineProperty(F, "__esModule", {
  value: !0
});
var L = F.default = void 0, ne = te(I()), re = G, ae = (0, ne.default)(/* @__PURE__ */ (0, re.jsx)("path", {
  d: "M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"
}), "GetApp");
L = F.default = ae;
const se = B()({
  dialogContent: {
    width: "80em"
  },
  textAreaFont: {
    fontFamily: "Courier New"
  }
});
async function oe(n, s, p) {
  const v = b.getSession(n), { leftOffset: a, rightOffset: c } = n;
  if (!a || !c)
    throw new Error("no offsets on model to use for range");
  if (a.assemblyName !== c.assemblyName)
    throw new Error("not able to fetch sequences from multiple assemblies");
  const { rpcManager: i, assemblyManager: E } = v, l = a.assemblyName || c.assemblyName || "", h = E.get(l);
  if (!h)
    throw new Error(`assembly ${l} not found`);
  const C = Z.getConf(h, ["sequence", "adapter"]), u = "getSequence";
  return i.call(u, "CoreGetFeatures", {
    adapterConfig: C,
    regions: s,
    sessionId: u,
    signal: p
  });
}
function le({ model: n, handleClose: s }) {
  const { classes: p } = se(), v = b.getSession(n), [a, c] = d(), [i, E] = d(), [l, h] = d(!1), [C, u] = d(!1), [q, O] = d(!1), { leftOffset: _, rightOffset: R } = n, y = i === void 0;
  T(() => {
    let t = !0;
    const r = new AbortController();
    return (async () => {
      try {
        const o = n.getSelectedRegions(_, R);
        if (o.length === 0)
          throw new Error("Selected region is out of bounds");
        const m = await oe(n, o, r.signal);
        t && E(m);
      } catch (o) {
        console.error(o), t && c(o);
      }
    })(), () => {
      r.abort(), t = !1;
    };
  }, [n, v, _, R]);
  const f = i ? j(i.filter((t) => !!t).map((t) => {
    let r = t.get("seq");
    const o = t.get("refName"), m = t.get("start") + 1, w = t.get("end"), k = `${o}:${m}-${w}`;
    if ((r == null ? void 0 : r.length) !== w - m + 1)
      throw new Error(`${k} returned ${r.length.toLocaleString()} bases, but should have returned ${(w - m).toLocaleString()}`);
    return l && (r = b.reverse(r)), q && (r = b.complement(r)), {
      header: k + (l ? "-rev" : "") + (q ? "-comp" : ""),
      seq: r
    };
  })) : "", S = f ? f.length > 1e6 : !1;
  return e.createElement(
    P.Dialog,
    { maxWidth: "xl", open: !0, onClose: () => {
      s(), n.setOffsets();
    }, title: "Reference sequence" },
    e.createElement(
      H,
      null,
      a ? e.createElement(A, { color: "error" }, `${a}`) : y ? e.createElement(
        W,
        null,
        "Retrieving reference sequence...",
        e.createElement(V, { style: { marginLeft: 10 }, size: 20, disableShrink: !0 })
      ) : null,
      e.createElement(J, { "data-testid": "rubberband-sequence", variant: "outlined", multiline: !0, minRows: 5, maxRows: 10, disabled: S, className: p.dialogContent, fullWidth: !0, value: S ? "Reference sequence too large to display, use the download FASTA button" : f, InputProps: {
        readOnly: !0,
        classes: {
          input: p.textAreaFont
        }
      } }),
      e.createElement(
        K,
        null,
        e.createElement($, { control: e.createElement(N, { value: l, onChange: (t) => h(t.target.checked) }), label: "Reverse sequence" }),
        e.createElement($, { control: e.createElement(N, { value: q, onChange: (t) => O(t.target.checked) }), label: "Complement sequence" })
      ),
      e.createElement(A, { style: { margin: 10 } }, 'Note: Check both boxes for the "reverse complement"')
    ),
    e.createElement(
      Q,
      null,
      e.createElement(x, { onClick: () => {
        U(f), u(!0), setTimeout(() => u(!1), 500);
      }, disabled: y || !!a || S, color: "primary", startIcon: e.createElement(X, null) }, C ? "Copied" : "Copy to clipboard"),
      e.createElement(x, { onClick: () => {
        Y.saveAs(new Blob([f || ""], {
          type: "text/x-fasta;charset=utf-8"
        }), "jbrowse_ref_seq.fa");
      }, disabled: y || !!a, color: "primary", startIcon: e.createElement(L, null) }, "Download FASTA"),
      e.createElement(x, { onClick: s, variant: "contained" }, "Close")
    )
  );
}
const fe = z(le);
export {
  fe as default
};
