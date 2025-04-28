import e, { useState as m } from "react";
import { o as u, q as h, u as g, b as d, T as p, n as x, c as y, B as n } from "./index-8ee3baa4.js";
import "react-dom";
const f = h()({
  root: {
    width: 500
  }
});
function E(l) {
  const { model: o, handleClose: a } = l, { classes: r } = f(), { maxHeight: s = "" } = o, [t, i] = m(`${s}`);
  return e.createElement(
    g.Dialog,
    { open: !0, onClose: a, title: "Filter options" },
    e.createElement(
      d,
      { className: r.root },
      e.createElement(p, null, 'Set max height for the track. For example, you can increase this if the layout says "Max height reached"'),
      e.createElement(x, { value: t, onChange: (c) => i(c.target.value), placeholder: "Enter max height for layout" }),
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
const v = u(E);
export {
  v as default
};
