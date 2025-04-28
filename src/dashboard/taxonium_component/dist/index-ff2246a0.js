import e from "react";
import { o as r, D as l, a as i, b as c, T as u, c as f, B as a } from "./index-8ee3baa4.js";
import "react-dom";
function m(t) {
  return typeof t == "object" && !!t && "jbrowse" in t;
}
const d = r(({ rootModel: t, onClose: n }) => {
  if (!t)
    return null;
  if (!m(t))
    return console.error("Incorrect rootmodel"), null;
  const { jbrowse: s, session: o } = t;
  return e.createElement(
    l,
    { open: !0, onClose: n },
    e.createElement(i, null, "Set default session"),
    e.createElement(
      c,
      null,
      e.createElement(u, null, 'Select "Set current session as default" to make your current session saved to the config file. You can also hit "Clear default session", which would remove the default session from the config.')
    ),
    e.createElement(
      f,
      null,
      e.createElement(a, { variant: "contained", onClick: () => {
        s.setDefaultSessionConf({
          name: "New session"
        }), n();
      } }, "Clear default session"),
      e.createElement(a, { color: "secondary", variant: "contained", onClick: () => n() }, "Cancel"),
      e.createElement(a, { color: "primary", variant: "contained", onClick: () => {
        s.setDefaultSessionConf(o), n();
      } }, "Set current session as default")
    )
  );
}), E = d;
export {
  E as default
};
