import t, { useState as E, useEffect as P, lazy as se } from "react";
import { a1 as ce, o as d, a2 as j, a3 as ue, a4 as R, n as v, a5 as k, I as $, m as F, B as A, a6 as N, q as w, t as D, a7 as me, a8 as de, a9 as pe, aa as he, ab as S, ac as T, ad as H, p as U, u as B, ae as ge, s as fe, C as Ee, af as ve, ag as _, ah as q, M as G, P as J, ai as $e, aj as be, ak as W, l as x, al as K, am as X, an as Z, T as Q, ao as Y, F as ye, ap as Ce } from "./index-8ee3baa4.js";
import "react-dom";
function xe(e) {
  throw new Error('Could not dynamically require "' + e + '". Please configure the dynamicRequireTargets or/and ignoreDynamicRequires option of @rollup/plugin-commonjs appropriately for this require call to work.');
}
var ee = { exports: {} };
(function(e, a) {
  (function(n, o) {
    typeof xe == "function" ? e.exports = o() : n.pluralize = o();
  })(ce, function() {
    var n = [], o = [], r = {}, s = {}, c = {};
    function u(l) {
      return typeof l == "string" ? new RegExp("^" + l + "$", "i") : l;
    }
    function p(l, i) {
      return l === i ? i : l === l.toLowerCase() ? i.toLowerCase() : l === l.toUpperCase() ? i.toUpperCase() : l[0] === l[0].toUpperCase() ? i.charAt(0).toUpperCase() + i.substr(1).toLowerCase() : i.toLowerCase();
    }
    function b(l, i) {
      return l.replace(/\$(\d{1,2})/g, function(f, h) {
        return i[h] || "";
      });
    }
    function I(l, i) {
      return l.replace(i[0], function(f, h) {
        var g = b(i[1], arguments);
        return p(f === "" ? l[h - 1] : f, g);
      });
    }
    function y(l, i, f) {
      if (!l.length || r.hasOwnProperty(l))
        return i;
      for (var h = f.length; h--; ) {
        var g = f[h];
        if (g[0].test(i))
          return I(i, g);
      }
      return i;
    }
    function M(l, i, f) {
      return function(h) {
        var g = h.toLowerCase();
        return i.hasOwnProperty(g) ? p(h, g) : l.hasOwnProperty(g) ? p(h, l[g]) : y(g, h, f);
      };
    }
    function z(l, i, f, h) {
      return function(g) {
        var C = g.toLowerCase();
        return i.hasOwnProperty(C) ? !0 : l.hasOwnProperty(C) ? !1 : y(C, C, f) === C;
      };
    }
    function m(l, i, f) {
      var h = i === 1 ? m.singular(l) : m.plural(l);
      return (f ? i + " " : "") + h;
    }
    return m.plural = M(
      c,
      s,
      n
    ), m.isPlural = z(
      c,
      s,
      n
    ), m.singular = M(
      s,
      c,
      o
    ), m.isSingular = z(
      s,
      c,
      o
    ), m.addPluralRule = function(l, i) {
      n.push([u(l), i]);
    }, m.addSingularRule = function(l, i) {
      o.push([u(l), i]);
    }, m.addUncountableRule = function(l) {
      if (typeof l == "string") {
        r[l.toLowerCase()] = !0;
        return;
      }
      m.addPluralRule(l, "$0"), m.addSingularRule(l, "$0");
    }, m.addIrregularRule = function(l, i) {
      i = i.toLowerCase(), l = l.toLowerCase(), c[l] = i, s[i] = l;
    }, [
      // Pronouns.
      ["I", "we"],
      ["me", "us"],
      ["he", "they"],
      ["she", "they"],
      ["them", "them"],
      ["myself", "ourselves"],
      ["yourself", "yourselves"],
      ["itself", "themselves"],
      ["herself", "themselves"],
      ["himself", "themselves"],
      ["themself", "themselves"],
      ["is", "are"],
      ["was", "were"],
      ["has", "have"],
      ["this", "these"],
      ["that", "those"],
      // Words ending in with a consonant and `o`.
      ["echo", "echoes"],
      ["dingo", "dingoes"],
      ["volcano", "volcanoes"],
      ["tornado", "tornadoes"],
      ["torpedo", "torpedoes"],
      // Ends with `us`.
      ["genus", "genera"],
      ["viscus", "viscera"],
      // Ends with `ma`.
      ["stigma", "stigmata"],
      ["stoma", "stomata"],
      ["dogma", "dogmata"],
      ["lemma", "lemmata"],
      ["schema", "schemata"],
      ["anathema", "anathemata"],
      // Other irregular rules.
      ["ox", "oxen"],
      ["axe", "axes"],
      ["die", "dice"],
      ["yes", "yeses"],
      ["foot", "feet"],
      ["eave", "eaves"],
      ["goose", "geese"],
      ["tooth", "teeth"],
      ["quiz", "quizzes"],
      ["human", "humans"],
      ["proof", "proofs"],
      ["carve", "carves"],
      ["valve", "valves"],
      ["looey", "looies"],
      ["thief", "thieves"],
      ["groove", "grooves"],
      ["pickaxe", "pickaxes"],
      ["passerby", "passersby"]
    ].forEach(function(l) {
      return m.addIrregularRule(l[0], l[1]);
    }), [
      [/s?$/i, "s"],
      [/[^\u0000-\u007F]$/i, "$0"],
      [/([^aeiou]ese)$/i, "$1"],
      [/(ax|test)is$/i, "$1es"],
      [/(alias|[^aou]us|t[lm]as|gas|ris)$/i, "$1es"],
      [/(e[mn]u)s?$/i, "$1s"],
      [/([^l]ias|[aeiou]las|[ejzr]as|[iu]am)$/i, "$1"],
      [/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i, "$1i"],
      [/(alumn|alg|vertebr)(?:a|ae)$/i, "$1ae"],
      [/(seraph|cherub)(?:im)?$/i, "$1im"],
      [/(her|at|gr)o$/i, "$1oes"],
      [/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|automat|quor)(?:a|um)$/i, "$1a"],
      [/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)(?:a|on)$/i, "$1a"],
      [/sis$/i, "ses"],
      [/(?:(kni|wi|li)fe|(ar|l|ea|eo|oa|hoo)f)$/i, "$1$2ves"],
      [/([^aeiouy]|qu)y$/i, "$1ies"],
      [/([^ch][ieo][ln])ey$/i, "$1ies"],
      [/(x|ch|ss|sh|zz)$/i, "$1es"],
      [/(matr|cod|mur|sil|vert|ind|append)(?:ix|ex)$/i, "$1ices"],
      [/\b((?:tit)?m|l)(?:ice|ouse)$/i, "$1ice"],
      [/(pe)(?:rson|ople)$/i, "$1ople"],
      [/(child)(?:ren)?$/i, "$1ren"],
      [/eaux$/i, "$0"],
      [/m[ae]n$/i, "men"],
      ["thou", "you"]
    ].forEach(function(l) {
      return m.addPluralRule(l[0], l[1]);
    }), [
      [/s$/i, ""],
      [/(ss)$/i, "$1"],
      [/(wi|kni|(?:after|half|high|low|mid|non|night|[^\w]|^)li)ves$/i, "$1fe"],
      [/(ar|(?:wo|[ae])l|[eo][ao])ves$/i, "$1f"],
      [/ies$/i, "y"],
      [/\b([pl]|zomb|(?:neck|cross)?t|coll|faer|food|gen|goon|group|lass|talk|goal|cut)ies$/i, "$1ie"],
      [/\b(mon|smil)ies$/i, "$1ey"],
      [/\b((?:tit)?m|l)ice$/i, "$1ouse"],
      [/(seraph|cherub)im$/i, "$1"],
      [/(x|ch|ss|sh|zz|tto|go|cho|alias|[^aou]us|t[lm]as|gas|(?:her|at|gr)o|[aeiou]ris)(?:es)?$/i, "$1"],
      [/(analy|diagno|parenthe|progno|synop|the|empha|cri|ne)(?:sis|ses)$/i, "$1sis"],
      [/(movie|twelve|abuse|e[mn]u)s$/i, "$1"],
      [/(test)(?:is|es)$/i, "$1is"],
      [/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i, "$1us"],
      [/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|quor)a$/i, "$1um"],
      [/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)a$/i, "$1on"],
      [/(alumn|alg|vertebr)ae$/i, "$1a"],
      [/(cod|mur|sil|vert|ind)ices$/i, "$1ex"],
      [/(matr|append)ices$/i, "$1ix"],
      [/(pe)(rson|ople)$/i, "$1rson"],
      [/(child)ren$/i, "$1"],
      [/(eau)x?$/i, "$1"],
      [/men$/i, "man"]
    ].forEach(function(l) {
      return m.addSingularRule(l[0], l[1]);
    }), [
      // Singular words with no plurals.
      "adulthood",
      "advice",
      "agenda",
      "aid",
      "aircraft",
      "alcohol",
      "ammo",
      "analytics",
      "anime",
      "athletics",
      "audio",
      "bison",
      "blood",
      "bream",
      "buffalo",
      "butter",
      "carp",
      "cash",
      "chassis",
      "chess",
      "clothing",
      "cod",
      "commerce",
      "cooperation",
      "corps",
      "debris",
      "diabetes",
      "digestion",
      "elk",
      "energy",
      "equipment",
      "excretion",
      "expertise",
      "firmware",
      "flounder",
      "fun",
      "gallows",
      "garbage",
      "graffiti",
      "hardware",
      "headquarters",
      "health",
      "herpes",
      "highjinks",
      "homework",
      "housework",
      "information",
      "jeans",
      "justice",
      "kudos",
      "labour",
      "literature",
      "machinery",
      "mackerel",
      "mail",
      "media",
      "mews",
      "moose",
      "music",
      "mud",
      "manga",
      "news",
      "only",
      "personnel",
      "pike",
      "plankton",
      "pliers",
      "police",
      "pollution",
      "premises",
      "rain",
      "research",
      "rice",
      "salmon",
      "scissors",
      "series",
      "sewage",
      "shambles",
      "shrimp",
      "software",
      "species",
      "staff",
      "swine",
      "tennis",
      "traffic",
      "transportation",
      "trout",
      "tuna",
      "wealth",
      "welfare",
      "whiting",
      "wildebeest",
      "wildlife",
      "you",
      /pok[eé]mon$/i,
      // Regexes.
      /[^aeiou]ese$/i,
      // "chinese", "japanese"
      /deer$/i,
      // "deer", "reindeer"
      /fish$/i,
      // "fish", "blowfish", "angelfish"
      /measles$/i,
      /o[iu]s$/i,
      // "carnivorous"
      /pox$/i,
      // "chickpox", "smallpox"
      /sheep$/i
    ].forEach(m.addUncountableRule), m;
  });
})(ee);
var we = ee.exports;
const te = d(({ slot: e }) => {
  const [a, n] = E(""), [o, r] = E(!1);
  return t.createElement(
    t.Fragment,
    null,
    e.name ? t.createElement(j, null, e.name) : null,
    t.createElement(
      ue,
      { disablePadding: !0 },
      e.value.map((s, c) => t.createElement(
        R,
        { key: c, disableGutters: !0 },
        t.createElement(v, { value: s, onChange: (u) => e.setAtIndex(c, u.target.value), InputProps: {
          endAdornment: t.createElement(
            k,
            { position: "end" },
            t.createElement(
              $,
              { onClick: () => e.removeAtIndex(c) },
              t.createElement(F, null)
            )
          )
        } })
      )),
      o ? t.createElement(
        R,
        { disableGutters: !0 },
        t.createElement(v, { value: a, placeholder: "add new", onChange: (s) => n(s.target.value), InputProps: {
          endAdornment: t.createElement(
            k,
            { position: "end" },
            t.createElement(
              t.Fragment,
              null,
              t.createElement(A, { color: "primary", variant: "contained", style: { margin: 2 }, "data-testid": `stringArrayAdd-${e.name}`, onClick: () => {
                r(!1), e.add(a), n("");
              } }, "OK"),
              t.createElement(A, { color: "primary", variant: "contained", style: { margin: 2 }, onClick: () => {
                r(!1), n("");
              } }, "Cancel")
            )
          )
        } })
      ) : null,
      t.createElement(A, { color: "primary", variant: "contained", style: { margin: 4 }, disabled: o, onClick: () => r(!0) }, "Add item")
    ),
    t.createElement(N, null, e.description)
  );
}), V = 'Consolas, "Andale Mono WT", "Andale Mono", "Lucida Console", "Lucida Sans Typewriter", "DejaVu Sans Mono", "Bitstream Vera Sans Mono", "Liberation Mono", "Nimbus Mono L", Monaco, "Courier New", Courier, monospace', ke = w()((e) => ({
  callbackEditor: {
    marginTop: "16px",
    borderBottom: `1px solid ${e.palette.divider}`,
    width: "100%",
    fontFamily: V
  },
  textAreaFont: {
    fontFamily: V
  },
  callbackContainer: {
    width: "100%",
    overflowX: "auto"
  },
  error: {
    color: "red",
    fontSize: "0.8em"
  }
}));
function Se({ slot: e }) {
  const { classes: a } = ke(), [n, o] = E(e.value), [r, s] = E(), c = D.useDebounce(n, 400);
  return P(() => {
    var u;
    try {
      const p = c.startsWith("jexl:") ? c : `jexl:${c}`;
      if (p === "jexl:")
        throw new Error("Empty jexl expression is not valid");
      me(p, (u = de(e).pluginManager) === null || u === void 0 ? void 0 : u.jexl), e.set(p), s(void 0);
    } catch (p) {
      console.error({ e: p }), s(p);
    }
  }, [c, e]), t.createElement(
    t.Fragment,
    null,
    r ? t.createElement("p", { className: a.error }, `${r}`) : null,
    t.createElement(
      "div",
      { className: a.callbackContainer },
      t.createElement(v, { multiline: !0, className: a.callbackEditor, value: n.startsWith("jexl:") ? n.split("jexl:")[1] : n, onChange: (u) => o(u.target.value), style: { background: r ? "#fdd" : void 0 }, InputProps: {
        classes: {
          input: a.textAreaFont
        }
      } }),
      t.createElement("p", null, e.description),
      t.createElement(
        pe,
        { title: t.createElement(
          "div",
          null,
          "Callbacks are written in Jexl format. Click to learn more.",
          t.createElement("br", null),
          " Names of available context items: ",
          e.contextVariable
        ), arrow: !0 },
        t.createElement(
          $,
          { color: "primary", onClick: () => {
            const u = window.open("https://github.com/TomFrost/Jexl", "_blank", "noopener,noreferrer");
            u && (u.opener = null);
          } },
          t.createElement(he, null)
        )
      )
    )
  );
}
const Te = d(Se), Ne = se(() => import("./ColorPicker-affa19c8.js").then((e) => e.C)), Le = (e) => {
  const { value: a = "#000", label: n = "", TextFieldProps: o = {}, onChange: r } = e, [s, c] = E(!1);
  return t.createElement(
    "div",
    { style: { display: "flex" } },
    t.createElement(v, { value: a, label: n, onClick: () => c(!s), onChange: (u) => r(u.target.value), ...o }),
    t.createElement(
      "div",
      { style: { marginTop: 10 } },
      t.createElement(
        t.Suspense,
        { fallback: t.createElement("div", null) },
        t.createElement(Ne, { color: a, onChange: (u) => r(u) })
      )
    )
  );
};
function Ie(e) {
  const { slot: a } = e;
  return t.createElement(Le, { label: a.name, value: a.value, onChange: (n) => a.set(n), TextFieldProps: {
    helperText: a.description,
    fullWidth: !0
  } });
}
const Ae = d(Ie), Pe = w()((e) => ({
  card: {
    marginTop: e.spacing(1)
  }
})), je = d(function({ slot: e }) {
  const { classes: a } = Pe(), [n, o] = E("");
  return t.createElement(
    t.Fragment,
    null,
    t.createElement(j, null, e.name),
    [...e.value].map(([r, s]) => t.createElement(
      S,
      { raised: !0, key: r, className: a.card },
      t.createElement(T, { title: r, action: t.createElement(
        $,
        { onClick: () => e.remove(r) },
        t.createElement(F, null)
      ) }),
      t.createElement(
        H,
        null,
        t.createElement(te, { slot: {
          name: e.name,
          value: s,
          description: `Values associated with entry ${r}`,
          setAtIndex: (c, u) => e.setAtKeyIndex(r, c, u),
          removeAtIndex: (c) => e.removeAtKeyIndex(r, c),
          add: (c) => e.addToKey(r, c)
        } })
      )
    )),
    t.createElement(
      S,
      { raised: !0, className: a.card },
      t.createElement(T, { disableTypography: !0, title: t.createElement(v, { fullWidth: !0, value: n, placeholder: "add new", onChange: (r) => o(r.target.value), InputProps: {
        endAdornment: t.createElement(
          k,
          { position: "end" },
          t.createElement(
            $,
            { disabled: n === "", onClick: () => {
              e.add(n, []), o("");
            } },
            t.createElement(U, null)
          )
        )
      } }) })
    ),
    t.createElement(N, null, e.description)
  );
});
function L(e) {
  const { helperText: a } = e;
  return t.createElement(v, { ...e, helperText: t.createElement(B.SanitizedHTML, { html: a || "" }), FormHelperTextProps: {
    // @ts-expect-error
    component: "div"
  }, fullWidth: !0 });
}
const ae = d(function({ slot: e }) {
  const [a, n] = E(e.value);
  return P(() => {
    var o;
    const r = Number.parseFloat(a);
    Number.isNaN(r) ? (o = e.reset) === null || o === void 0 || o.call(e) : e.set(r);
  }, [e, a]), t.createElement(L, { label: e.name, helperText: e.description, value: a, type: "number", onChange: (o) => n(o.target.value) });
}), Fe = w()((e) => ({
  card: {
    marginTop: e.spacing(1)
  }
})), Me = d(function({ slot: e }) {
  const { classes: a } = Fe(), [n, o] = E("");
  return t.createElement(
    t.Fragment,
    null,
    t.createElement(j, null, e.name),
    [...e.value].map(([r, s]) => t.createElement(
      S,
      { raised: !0, key: r, className: a.card },
      t.createElement(T, { title: r, action: t.createElement(
        $,
        { onClick: () => e.remove(r) },
        t.createElement(F, null)
      ) }),
      t.createElement(
        H,
        null,
        t.createElement(ae, { slot: {
          value: s,
          set: (c) => e.add(r, c)
        } })
      )
    )),
    t.createElement(
      S,
      { raised: !0, className: a.card },
      t.createElement(T, { disableTypography: !0, title: t.createElement(v, { fullWidth: !0, value: n, placeholder: "add new", onChange: (r) => o(r.target.value), InputProps: {
        endAdornment: t.createElement(
          k,
          { position: "end" },
          t.createElement(
            $,
            { disabled: n === "", onClick: () => {
              e.add(n, 0), o("");
            } },
            t.createElement(U, null)
          )
        )
      } }) })
    ),
    t.createElement(N, null, e.description)
  );
}), ze = d(function({ slot: e }) {
  return t.createElement(
    ge,
    null,
    t.createElement(fe, { label: e.name, control: t.createElement(Ee, { checked: e.value, onChange: (a) => e.set(a.target.checked) }) }),
    t.createElement(N, null, e.description)
  );
}), ne = d(({ slot: e }) => t.createElement(L, { label: e.name, helperText: e.description, value: e.value, onChange: (a) => e.set(a.target.value) })), Re = d(({ slot: e }) => t.createElement(v, { label: e.name, helperText: e.description, multiline: !0, value: e.value, onChange: (a) => e.set(a.target.value) })), _e = () => t.createElement(
  be,
  null,
  t.createElement("path", { d: "M20.41,3C21.8,5.71 22.35,8.84 22,12C21.8,15.16 20.7,18.29 18.83,21L17.3,20C18.91,17.57 19.85,14.8 20,12C20.34,9.2 19.89,6.43 18.7,4L20.41,3M5.17,3L6.7,4C5.09,6.43 4.15,9.2 4,12C3.66,14.8 4.12,17.57 5.3,20L3.61,21C2.21,18.29 1.65,15.17 2,12C2.2,8.84 3.3,5.71 5.17,3M12.08,10.68L14.4,7.45H16.93L13.15,12.45L15.35,17.37H13.09L11.71,14L9.28,17.33H6.76L10.66,12.21L8.53,7.45H10.8L12.08,10.68Z" })
), qe = d(({ slot: e }) => {
  const [a, n] = E(e.value);
  return P(() => {
    const o = Number.parseInt(a, 10);
    Number.isNaN(o) || e.set(o);
  }, [e, a]), t.createElement(L, { label: e.name, helperText: e.description, value: a, type: "number", onChange: (o) => n(o.target.value) });
}), We = d(function({ slot: e, slotSchema: a }) {
  const n = ve(_(a)), o = q(q(_(n.properties.value))[1]).map((r) => r.value);
  return t.createElement(L, { value: e.value, label: e.name, select: !0, helperText: e.description, onChange: (r) => e.set(r.target.value) }, o.map((r) => t.createElement(G, { key: r, value: r }, r)));
}), Ve = d(function({ slot: e }) {
  var a;
  return t.createElement(B.FileSelector, { location: e.value, setLocation: (n) => e.set(n), name: e.name, description: e.description, rootModel: (a = D.getEnv(e).pluginManager) === null || a === void 0 ? void 0 : a.rootModel });
}), O = {
  string: ne,
  text: Re,
  fileLocation: Ve,
  stringArray: te,
  stringArrayMap: je,
  numberMap: Me,
  number: ae,
  integer: qe,
  color: Ae,
  stringEnum: We,
  boolean: ze,
  frozen: W,
  configRelationships: W
}, re = w()((e) => ({
  paper: {
    display: "flex",
    marginBottom: e.spacing(2),
    position: "relative"
  },
  paperContent: {
    width: "100%"
  },
  slotModeSwitch: {
    width: 24,
    background: e.palette.secondary.light,
    display: "flex",
    justifyContent: "center",
    alignItems: "center"
  }
})), Oe = d(
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  ({ slot: e, slotSchema: a }) => {
    const { classes: n } = re(), { type: o } = e;
    let r = e.isCallback ? Te : (
      // @ts-expect-error
      O[o]
    );
    return r || (console.warn(`no slot editor defined for ${o}, editing as string`), r = ne), o in O || console.warn(`SlotEditor needs to implement ${o}`), t.createElement(
      J,
      { className: n.paper },
      t.createElement(
        "div",
        { className: n.paperContent },
        t.createElement(r, { slot: e, slotSchema: a })
      ),
      t.createElement("div", { className: n.slotModeSwitch }, e.contextVariable.length ? t.createElement($, { onClick: () => e.isCallback ? e.convertToValue() : e.convertToCallback(), title: `convert to ${e.isCallback ? "regular value" : "callback"}` }, e.isCallback ? t.createElement(_e, null) : t.createElement($e, null)) : null)
    );
  }
), De = d(({ typeNameChoices: e, slot: a, slotName: n, onChange: o }) => {
  const { classes: r } = re();
  return t.createElement(
    J,
    { className: r.paper },
    t.createElement(
      "div",
      { className: r.paperContent },
      t.createElement(v, { value: a.type, label: "Type", select: !0, helperText: `Type of ${n} to use`, fullWidth: !0, onChange: o }, e.map((s) => t.createElement(G, { key: s, value: s }, s)))
    )
  );
}), le = w()((e) => {
  var a;
  return {
    expandIcon: {
      color: ((a = e.palette.tertiary) === null || a === void 0 ? void 0 : a.contrastText) || "#fff"
    },
    root: {
      padding: e.spacing(1, 3, 1, 1)
    },
    expansionPanelDetails: {
      display: "block",
      padding: e.spacing(1)
    },
    accordion: {
      border: `1px solid ${e.palette.text.primary}`
    },
    noOverflow: {
      width: "100%",
      overflowX: "auto"
    }
  };
}), oe = d(function(e) {
  const { classes: a } = le(), { slotName: n, slotSchema: o, schema: r, slot: s = r[n], path: c = [] } = e;
  let u;
  if (x.isConfigurationSchemaType(o)) {
    if (s.length)
      return s.map((b, I) => {
        const y = `${we.singular(n)} ${I + 1}`;
        return t.createElement(oe, { ...e, key: y, slot: b, slotName: y });
      });
    const p = x.getTypeNamesFromExplicitlyTypedUnion(o);
    return p.length && (u = t.createElement(De, { typeNameChoices: p, slotName: n, slot: s, onChange: (b) => {
      b.target.value !== s.type && r.setSubschema(n, { type: b.target.value });
    } })), t.createElement(
      K,
      { defaultExpanded: !0, className: a.accordion },
      t.createElement(
        X,
        { expandIcon: t.createElement(Z, { className: a.expandIcon }) },
        t.createElement(Q, null, [...c, n].join("➔"))
      ),
      t.createElement(
        Y,
        { className: a.expansionPanelDetails },
        u,
        t.createElement(
          ye,
          { className: a.noOverflow },
          t.createElement(ie, { schema: s, path: [...c, n] })
        )
      )
    );
  }
  return x.isConfigurationSlotType(o) ? t.createElement(Oe, { key: n, slot: s, slotSchema: o }) : null;
}), ie = d(function({ schema: e, path: a = [] }) {
  const n = Ce(e).properties;
  return t.createElement(t.Fragment, null, Object.entries(n).map(([o, r]) => t.createElement(oe, { key: o, slotName: o, slotSchema: r, path: a, schema: e })));
}), Ge = d(function({ model: e }) {
  const { classes: a } = le(), { target: n } = e, o = n && x.readConfObject(n, "trackId"), r = n && x.readConfObject(n, "name");
  return t.createElement(
    K,
    { key: o, defaultExpanded: !0, className: a.accordion },
    t.createElement(
      X,
      { expandIcon: t.createElement(Z, { className: a.expandIcon }) },
      t.createElement(Q, null, r ?? "Configuration")
    ),
    t.createElement(Y, { className: a.expansionPanelDetails, "data-testid": "configEditor" }, n ? t.createElement(ie, { schema: n }) : "no target set")
  );
});
export {
  Ge as default
};
