import e, { useState as C } from "react";
import { o as m, u as g, b as E, T as l, s as i, R as s, c as p, B as u } from "./index-8ee3baa4.js";
import { a as n } from "./ColorPicker-affa19c8.js";
import "react-dom";
import "./index-090e001f.js";
function f({ model: o, handleClose: a }) {
  const [r, c] = C(!1);
  return e.createElement(
    g.Dialog,
    { open: !0, onClose: a, title: "Set color" },
    e.createElement(
      E,
      null,
      e.createElement(l, null, "Select either an overall color, or the positive/negative colors. Note that density renderers only work properly with positive/negative colors"),
      e.createElement(i, { checked: !r, onClick: () => c(!1), control: e.createElement(s, null), label: "Overall color" }),
      e.createElement(i, { checked: r, onClick: () => c(!0), control: e.createElement(s, null), label: "Positive/negative color" }),
      r ? e.createElement(
        e.Fragment,
        null,
        e.createElement(l, null, "Positive color"),
        e.createElement(n.ColorPicker, { color: o.posColor || "black", onChange: (t) => {
          o.setPosColor(t), o.setColor(void 0);
        } }),
        e.createElement(l, null, "Negative color"),
        e.createElement(n.ColorPicker, { color: o.negColor || "black", onChange: (t) => {
          o.setNegColor(t), o.setColor(void 0);
        } })
      ) : e.createElement(
        e.Fragment,
        null,
        e.createElement(l, null, "Overall color"),
        e.createElement(n.ColorPicker, { color: o.color || "black", onChange: (t) => o.setColor(t) })
      )
    ),
    e.createElement(
      p,
      null,
      e.createElement(u, { onClick: () => {
        o.setPosColor(void 0), o.setNegColor(void 0), o.setColor(void 0);
      }, color: "secondary", variant: "contained" }, "Restore default"),
      e.createElement(u, { variant: "contained", color: "primary", type: "submit", onClick: () => a() }, "Submit")
    )
  );
}
const P = m(f);
export {
  P as default
};
