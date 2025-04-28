import e, { useState as l } from "react";
import { o as m, q as u, u as h, b as g, T as x, n as p, c as y, B as n } from "./index-8ee3baa4.js";
import "react-dom";
const d = u()({
  root: {
    width: 500
  }
});
function f({ model: o, handleClose: a }) {
  const { classes: r } = d(), { maxHeight: s = "" } = o, [t, c] = l(`${s}`);
  return e.createElement(
    h.Dialog,
    { open: !0, onClose: a, title: "Set max height" },
    e.createElement(
      g,
      { className: r.root },
      e.createElement(x, null, 'Set max height for the track. For example, you can increase this if the layout says "Max height reached"'),
      e.createElement(p, { value: t, onChange: (i) => c(i.target.value), placeholder: "Enter max score" }),
      e.createElement(
        y,
        null,
        e.createElement(n, { variant: "contained", color: "primary", type: "submit", autoFocus: !0, onClick: () => {
          o.setMaxHeight(t !== "" && !Number.isNaN(+t) ? +t : void 0), a();
        } }, "Submit"),
        e.createElement(n, { variant: "contained", color: "secondary", onClick: () => a() }, "Cancel")
      )
    )
  );
}
const v = m(f);
export {
  v as default
};
