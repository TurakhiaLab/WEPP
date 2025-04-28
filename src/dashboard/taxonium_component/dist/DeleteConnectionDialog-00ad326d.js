import e from "react";
import { o as r, D as m, a as u, b as D, aE as s, c as C, B as l } from "./index-8ee3baa4.js";
import "react-dom";
function E({ deleteDialogDetails: c, session: n, handleClose: o }) {
  const { connectionConf: a, name: i } = c;
  return e.createElement(
    m,
    { open: !0 },
    e.createElement(
      u,
      null,
      'Delete connection "',
      i,
      '"'
    ),
    e.createElement(
      D,
      null,
      e.createElement(s, null, "Are you sure you want to delete this connection?")
    ),
    e.createElement(
      C,
      null,
      e.createElement(l, { onClick: () => o(), color: "primary" }, "Cancel"),
      e.createElement(l, { variant: "contained", color: "primary", onClick: () => {
        var t;
        a && ((t = n.deleteConnection) === null || t === void 0 || t.call(n, a)), o();
      } }, "OK")
    )
  );
}
const d = r(E);
export {
  d as default
};
