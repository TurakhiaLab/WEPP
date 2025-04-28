import t, { useState as m } from "react";
import { u as p, b as x, T as o, n as u, c as f, B as M } from "./index-8ee3baa4.js";
import "react-dom";
function C(N) {
  const { model: a, handleClose: l } = N, { minScore: c, maxScore: s, scaleType: d } = a, [e, g] = m(`${c !== Number.MIN_VALUE ? c : ""}`), [r, E] = m(`${s !== Number.MAX_VALUE ? s : ""}`), i = e !== "" && r !== "" && !Number.isNaN(+e) && !Number.isNaN(+r) ? +r > +e : !0, b = d === "log" && e !== "" && !Number.isNaN(+e) ? +e > 0 : !0;
  return t.createElement(
    p.Dialog,
    { open: !0, onClose: l, title: "Set min/max score for track" },
    t.createElement(
      x,
      null,
      t.createElement(o, null, "Enter min/max score: "),
      i ? null : t.createElement(o, { color: "error" }, "Max is greater than or equal to min"),
      b ? null : t.createElement(o, { color: "error" }, "Min score should be greater than 0 for log scale"),
      t.createElement(u, { value: e, onChange: (n) => g(n.target.value), placeholder: "Enter min score" }),
      t.createElement(u, { value: r, onChange: (n) => E(n.target.value), placeholder: "Enter max score" })
    ),
    t.createElement(
      f,
      null,
      t.createElement(M, { variant: "contained", color: "primary", type: "submit", style: { marginLeft: 20 }, disabled: !i, onClick: () => {
        a.setMinScore(e !== "" && !Number.isNaN(+e) ? +e : void 0), a.setMaxScore(r !== "" && !Number.isNaN(+r) ? +r : void 0), l();
      } }, "Submit")
    )
  );
}
export {
  C as default
};
