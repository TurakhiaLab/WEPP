const f = class n {
  /**
   * Accept two comparable values and creates new instance of interval
   * Predicate Interval.comparable_less(low, high) supposed to return true on these values
   * @param low
   * @param high
   */
  constructor(e, t) {
    this.low = e, this.high = t;
  }
  /**
   * Clone interval
   * @returns {Interval}
   */
  clone() {
    return new n(this.low, this.high);
  }
  /**
   * Propery max returns clone of this interval
   * @returns {Interval}
   */
  get max() {
    return this.clone();
  }
  /**
   * Predicate returns true is this interval less than other interval
   * @param other_interval
   * @returns {boolean}
   */
  less_than(e) {
    return this.low < e.low || this.low == e.low && this.high < e.high;
  }
  /**
   * Predicate returns true is this interval equals to other interval
   * @param other_interval
   * @returns {boolean}
   */
  equal_to(e) {
    return this.low == e.low && this.high == e.high;
  }
  /**
   * Predicate returns true if this interval intersects other interval
   * @param other_interval
   * @returns {boolean}
   */
  intersect(e) {
    return !this.not_intersect(e);
  }
  /**
   * Predicate returns true if this interval does not intersect other interval
   * @param other_interval
   * @returns {boolean}
   */
  not_intersect(e) {
    return this.high < e.low || e.high < this.low;
  }
  /**
   * Returns new interval merged with other interval
   * @param {Interval} interval - Other interval to merge with
   * @returns {Interval}
   */
  merge(e) {
    return new n(
      this.low === void 0 ? e.low : Math.min(this.low, e.low),
      this.high === void 0 ? e.high : Math.max(this.high, e.high)
    );
  }
  /**
   * Returns how key should return
   */
  output() {
    return [this.low, this.high];
  }
  /**
   * Function returns maximum between two comparable values
   * @param interval1
   * @param interval2
   * @returns {Interval}
   */
  static comparable_max(e, t) {
    return e.merge(t);
  }
  /**
   * Predicate returns true if first value less than second value
   * @param val1
   * @param val2
   * @returns {boolean}
   */
  static comparable_less_than(e, t) {
    return e < t;
  }
}, s = 0, r = 1;
class h {
  constructor(e = void 0, t = void 0, i = null, l = null, a = null, u = r) {
    this.left = i, this.right = l, this.parent = a, this.color = u, this.item = { key: e, value: t }, e && e instanceof Array && e.length == 2 && !Number.isNaN(e[0]) && !Number.isNaN(e[1]) && (this.item.key = new f(Math.min(e[0], e[1]), Math.max(e[0], e[1]))), this.max = this.item.key ? this.item.key.max : void 0;
  }
  isNil() {
    return this.item.key === void 0 && this.item.value === void 0 && this.left === null && this.right === null && this.color === r;
  }
  less_than(e) {
    if (this.item.value === this.item.key && e.item.value === e.item.key)
      return this.item.key.less_than(e.item.key);
    {
      let t = this.item.value && e.item.value && this.item.value.less_than ? this.item.value.less_than(e.item.value) : this.item.value < e.item.value;
      return this.item.key.less_than(e.item.key) || this.item.key.equal_to(e.item.key) && t;
    }
  }
  equal_to(e) {
    if (this.item.value === this.item.key && e.item.value === e.item.key)
      return this.item.key.equal_to(e.item.key);
    {
      let t = this.item.value && e.item.value && this.item.value.equal_to ? this.item.value.equal_to(e.item.value) : this.item.value == e.item.value;
      return this.item.key.equal_to(e.item.key) && t;
    }
  }
  intersect(e) {
    return this.item.key.intersect(e.item.key);
  }
  copy_data(e) {
    this.item.key = e.item.key, this.item.value = e.item.value;
  }
  update_max() {
    if (this.max = this.item.key ? this.item.key.max : void 0, this.right && this.right.max) {
      const e = this.item.key.constructor.comparable_max;
      this.max = e(this.max, this.right.max);
    }
    if (this.left && this.left.max) {
      const e = this.item.key.constructor.comparable_max;
      this.max = e(this.max, this.left.max);
    }
  }
  // Other_node does not intersect any node of left subtree, if this.left.max < other_node.item.key.low
  not_intersect_left_subtree(e) {
    const t = this.item.key.constructor.comparable_less_than;
    let i = this.left.max.high !== void 0 ? this.left.max.high : this.left.max;
    return t(i, e.item.key.low);
  }
  // Other_node does not intersect right subtree if other_node.item.key.high < this.right.key.low
  not_intersect_right_subtree(e) {
    const t = this.item.key.constructor.comparable_less_than;
    let i = this.right.max.low !== void 0 ? this.right.max.low : this.right.item.key.low;
    return t(e.item.key.high, i);
  }
}
class o {
  /**
   * Construct new empty instance of IntervalTree
   */
  constructor() {
    this.root = null, this.nil_node = new h();
  }
  /**
   * Returns number of items stored in the interval tree
   * @returns {number}
   */
  get size() {
    let e = 0;
    return this.tree_walk(this.root, () => e++), e;
  }
  /**
   * Returns array of sorted keys in the ascending order
   * @returns {Array}
   */
  get keys() {
    let e = [];
    return this.tree_walk(this.root, (t) => e.push(
      t.item.key.output ? t.item.key.output() : t.item.key
    )), e;
  }
  /**
   * Return array of values in the ascending keys order
   * @returns {Array}
   */
  get values() {
    let e = [];
    return this.tree_walk(this.root, (t) => e.push(t.item.value)), e;
  }
  /**
   * Returns array of items (<key,value> pairs) in the ascended keys order
   * @returns {Array}
   */
  get items() {
    let e = [];
    return this.tree_walk(this.root, (t) => e.push({
      key: t.item.key.output ? t.item.key.output() : t.item.key,
      value: t.item.value
    })), e;
  }
  /**
   * Returns true if tree is empty
   * @returns {boolean}
   */
  isEmpty() {
    return this.root == null || this.root == this.nil_node;
  }
  /**
   * Clear tree
   */
  clear() {
    this.root = null;
  }
  /**
   * Insert new item into interval tree
   * @param {Interval} key - interval object or array of two numbers [low, high]
   * @param {any} value - value representing any object (optional)
   * @returns {Node} returns reference to inserted node as an object {key:interval, value: value}
   */
  insert(e, t = e) {
    if (e === void 0)
      return;
    let i = new h(e, t, this.nil_node, this.nil_node, null, s);
    return this.tree_insert(i), this.recalc_max(i), i;
  }
  /**
   * Returns true if item {key,value} exist in the tree
   * @param {Interval} key - interval correspondent to keys stored in the tree
   * @param {any} value - value object to be checked
   * @returns {boolean} true if item {key, value} exist in the tree, false otherwise
   */
  exist(e, t = e) {
    let i = new h(e, t);
    return !!this.tree_search(this.root, i);
  }
  /**
   * Remove entry {key, value} from the tree
   * @param {Interval} key - interval correspondent to keys stored in the tree
   * @param {any} value - value object
   * @returns {boolean} true if item {key, value} deleted, false if not found
   */
  remove(e, t = e) {
    let i = new h(e, t), l = this.tree_search(this.root, i);
    return l && this.tree_delete(l), l;
  }
  /**
   * Returns array of entry values which keys intersect with given interval <br/>
   * If no values stored in the tree, returns array of keys which intersect given interval
   * @param {Interval} interval - search interval, or tuple [low, high]
   * @param outputMapperFn(value,key) - optional function that maps (value, key) to custom output
   * @returns {Array}
   */
  search(e, t = (i, l) => i === l ? l.output() : i) {
    let i = new h(e), l = [];
    return this.tree_search_interval(this.root, i, l), l.map((a) => t(a.item.value, a.item.key));
  }
  /**
   * Returns true if intersection between given and any interval stored in the tree found
   * @param {Interval} interval - search interval or tuple [low, high]
   * @returns {boolean}
   */
  intersect_any(e) {
    let t = new h(e);
    return this.tree_find_any_interval(this.root, t);
  }
  /**
   * Tree visitor. For each node implement a callback function. <br/>
   * Method calls a callback function with two parameters (key, value)
   * @param visitor(key,value) - function to be called for each tree item
   */
  forEach(e) {
    this.tree_walk(this.root, (t) => e(t.item.key, t.item.value));
  }
  /** Value Mapper. Walk through every node and map node value to another value
  * @param callback(value,key) - function to be called for each tree item
  */
  map(e) {
    const t = new o();
    return this.tree_walk(this.root, (i) => t.insert(i.item.key, e(i.item.value, i.item.key))), t;
  }
  recalc_max(e) {
    let t = e;
    for (; t.parent != null; )
      t.parent.update_max(), t = t.parent;
  }
  tree_insert(e) {
    let t = this.root, i = null;
    if (this.root == null || this.root == this.nil_node)
      this.root = e;
    else {
      for (; t != this.nil_node; )
        i = t, e.less_than(t) ? t = t.left : t = t.right;
      e.parent = i, e.less_than(i) ? i.left = e : i.right = e;
    }
    this.insert_fixup(e);
  }
  // After insertion insert_node may have red-colored parent, and this is a single possible violation
  // Go upwords to the root and re-color until violation will be resolved
  insert_fixup(e) {
    let t, i;
    for (t = e; t != this.root && t.parent.color == s; )
      t.parent == t.parent.parent.left ? (i = t.parent.parent.right, i.color == s ? (t.parent.color = r, i.color = r, t.parent.parent.color = s, t = t.parent.parent) : (t == t.parent.right && (t = t.parent, this.rotate_left(t)), t.parent.color = r, t.parent.parent.color = s, this.rotate_right(t.parent.parent))) : (i = t.parent.parent.left, i.color == s ? (t.parent.color = r, i.color = r, t.parent.parent.color = s, t = t.parent.parent) : (t == t.parent.left && (t = t.parent, this.rotate_right(t)), t.parent.color = r, t.parent.parent.color = s, this.rotate_left(t.parent.parent)));
    this.root.color = r;
  }
  tree_delete(e) {
    let t, i;
    e.left == this.nil_node || e.right == this.nil_node ? t = e : t = this.tree_successor(e), t.left != this.nil_node ? i = t.left : i = t.right, i.parent = t.parent, t == this.root ? this.root = i : (t == t.parent.left ? t.parent.left = i : t.parent.right = i, t.parent.update_max()), this.recalc_max(i), t != e && (e.copy_data(t), e.update_max(), this.recalc_max(e)), /*fix_node != this.nil_node && */
    t.color == r && this.delete_fixup(i);
  }
  delete_fixup(e) {
    let t = e, i;
    for (; t != this.root && t.parent != null && t.color == r; )
      t == t.parent.left ? (i = t.parent.right, i.color == s && (i.color = r, t.parent.color = s, this.rotate_left(t.parent), i = t.parent.right), i.left.color == r && i.right.color == r ? (i.color = s, t = t.parent) : (i.right.color == r && (i.color = s, i.left.color = r, this.rotate_right(i), i = t.parent.right), i.color = t.parent.color, t.parent.color = r, i.right.color = r, this.rotate_left(t.parent), t = this.root)) : (i = t.parent.left, i.color == s && (i.color = r, t.parent.color = s, this.rotate_right(t.parent), i = t.parent.left), i.left.color == r && i.right.color == r ? (i.color = s, t = t.parent) : (i.left.color == r && (i.color = s, i.right.color = r, this.rotate_left(i), i = t.parent.left), i.color = t.parent.color, t.parent.color = r, i.left.color = r, this.rotate_right(t.parent), t = this.root));
    t.color = r;
  }
  tree_search(e, t) {
    if (!(e == null || e == this.nil_node))
      return t.equal_to(e) ? e : t.less_than(e) ? this.tree_search(e.left, t) : this.tree_search(e.right, t);
  }
  // Original search_interval method; container res support push() insertion
  // Search all intervals intersecting given one
  tree_search_interval(e, t, i) {
    e != null && e != this.nil_node && (e.left != this.nil_node && !e.not_intersect_left_subtree(t) && this.tree_search_interval(e.left, t, i), e.intersect(t) && i.push(e), e.right != this.nil_node && !e.not_intersect_right_subtree(t) && this.tree_search_interval(e.right, t, i));
  }
  tree_find_any_interval(e, t) {
    let i = !1;
    return e != null && e != this.nil_node && (e.left != this.nil_node && !e.not_intersect_left_subtree(t) && (i = this.tree_find_any_interval(e.left, t)), i || (i = e.intersect(t)), !i && e.right != this.nil_node && !e.not_intersect_right_subtree(t) && (i = this.tree_find_any_interval(e.right, t))), i;
  }
  local_minimum(e) {
    let t = e;
    for (; t.left != null && t.left != this.nil_node; )
      t = t.left;
    return t;
  }
  // not in use
  local_maximum(e) {
    let t = e;
    for (; t.right != null && t.right != this.nil_node; )
      t = t.right;
    return t;
  }
  tree_successor(e) {
    let t, i, l;
    if (e.right != this.nil_node)
      t = this.local_minimum(e.right);
    else {
      for (i = e, l = e.parent; l != null && l.right == i; )
        i = l, l = l.parent;
      t = l;
    }
    return t;
  }
  //           |            right-rotate(T,y)       |
  //           y            ---------------.       x
  //          / \                                  / \
  //         x   c          left-rotate(T,x)      a   y
  //        / \             <---------------         / \
  //       a   b                                    b   c
  rotate_left(e) {
    let t = e.right;
    e.right = t.left, t.left != this.nil_node && (t.left.parent = e), t.parent = e.parent, e == this.root ? this.root = t : e == e.parent.left ? e.parent.left = t : e.parent.right = t, t.left = e, e.parent = t, e != null && e != this.nil_node && e.update_max(), t = e.parent, t != null && t != this.nil_node && t.update_max();
  }
  rotate_right(e) {
    let t = e.left;
    e.left = t.right, t.right != this.nil_node && (t.right.parent = e), t.parent = e.parent, e == this.root ? this.root = t : e == e.parent.left ? e.parent.left = t : e.parent.right = t, t.right = e, e.parent = t, e != null && e != this.nil_node && e.update_max(), t = e.parent, t != null && t != this.nil_node && t.update_max();
  }
  tree_walk(e, t) {
    e != null && e != this.nil_node && (this.tree_walk(e.left, t), t(e), this.tree_walk(e.right, t));
  }
  /* Return true if all red nodes have exactly two black child nodes */
  testRedBlackProperty() {
    let e = !0;
    return this.tree_walk(this.root, function(t) {
      t.color == s && (t.left.color == r && t.right.color == r || (e = !1));
    }), e;
  }
  /* Throw error if not every path from root to bottom has same black height */
  testBlackHeightProperty(e) {
    let t = 0, i = 0, l = 0;
    if (e.color == r && t++, e.left != this.nil_node ? i = this.testBlackHeightProperty(e.left) : i = 1, e.right != this.nil_node ? l = this.testBlackHeightProperty(e.right) : l = 1, i != l)
      throw new Error("Red-black height property violated");
    return t += i, t;
  }
}
const m = o;
export {
  m as I
};
