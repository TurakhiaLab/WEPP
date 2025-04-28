import e, { useState as r } from "react";
import { o as p, u as d, b as S, T as b, n as f, s as C, C as y, c as E, B as c } from "./index-8ee3baa4.js";
import "react-dom";
function v(l) {
  const { model: a, handleClose: n } = l, { featureHeightSetting: s, noSpacing: u } = a, [t, g] = r(`${s}`), [i, h] = r(u), m = t !== "" && !Number.isNaN(+t);
  return e.createElement(
    d.Dialog,
    { open: !0, onClose: n, title: "Set feature height" },
    e.createElement(
      S,
      null,
      e.createElement(b, null, "Adjust the feature height and whether there is any spacing between features. Setting feature height to 1 and removing spacing makes the display very compact."),
      e.createElement(f, { value: t, helperText: "Feature height", onChange: (o) => g(o.target.value) }),
      e.createElement(C, { control: e.createElement(y, { checked: !!i, onChange: () => h((o) => !o) }), label: "Remove spacing between features in y-direction?" }),
      e.createElement(
        E,
        null,
        e.createElement(c, { variant: "contained", color: "primary", type: "submit", autoFocus: !0, disabled: !m, onClick: () => {
          a.setFeatureHeight(t !== "" && !Number.isNaN(+t) ? +t : void 0), a.setNoSpacing(i), n();
        } }, "Submit"),
        e.createElement(c, { variant: "contained", color: "secondary", onClick: () => n() }, "Cancel")
      )
    )
  );
}
const H = p(v);
export {
  H as default
};
