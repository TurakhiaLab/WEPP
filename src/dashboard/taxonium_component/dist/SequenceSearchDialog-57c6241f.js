import e, { useState as r } from "react";
import { o as q, q as f, u as x, b as T, T as p, n as w, F, s as o, C as l, c as A, B as S, t as D, v as _ } from "./index-8ee3baa4.js";
import "react-dom";
const N = f()({
  dialogContent: {
    width: "40em"
  }
});
function R({ model: n, handleClose: s }) {
  const { classes: C } = N(), [a, E] = r(""), [i, v] = r(!0), [u, y] = r(!0), [m, k] = r(!0);
  let c;
  try {
    new RegExp(a);
  } catch (t) {
    c = t;
  }
  return e.createElement(
    x.Dialog,
    { maxWidth: "xl", open: !0, onClose: s, title: "Sequence search" },
    e.createElement(
      T,
      { className: C.dialogContent },
      e.createElement(p, null, "Supply a sequence to search for. A track will be created with the resulting matches once submitted. You can also supply regex style expressions e.g. AACT(C|T)."),
      e.createElement(w, { value: a, onChange: (t) => E(t.target.value), helperText: "Sequence search pattern" }),
      e.createElement(
        F,
        null,
        e.createElement(o, { control: e.createElement(l, { checked: i, onChange: (t) => v(t.target.checked) }), label: "Search forward strand" }),
        e.createElement(o, { control: e.createElement(l, { checked: u, onChange: (t) => y(t.target.checked) }), label: "Search reverse strand" }),
        e.createElement(o, { control: e.createElement(l, { checked: m, onChange: (t) => k(t.target.checked) }), label: "Case insensitive" })
      ),
      c ? e.createElement(p, { color: "error" }, `${c}`) : null
    ),
    e.createElement(
      A,
      null,
      e.createElement(S, { onClick: () => {
        var t;
        if (a) {
          const h = `sequence_search_${+Date.now()}`, d = D.getSession(n), { assemblyManager: b } = d, g = n.assemblyNames[0];
          d.addTrackConf({
            trackId: h,
            name: `Sequence search ${a}`,
            assemblyNames: [g],
            type: "FeatureTrack",
            adapter: {
              type: "SequenceSearchAdapter",
              search: a,
              searchForward: i,
              searchReverse: u,
              caseInsensitive: m,
              sequenceAdapter: _((t = b.get(g)) === null || t === void 0 ? void 0 : t.configuration.sequence.adapter)
            }
          }), n.toggleTrack(h);
        }
        s();
      }, variant: "contained", color: "primary" }, "Submit"),
      e.createElement(S, { onClick: () => s(), variant: "contained", color: "secondary" }, "Close")
    )
  );
}
const G = q(R);
export {
  G as default
};
