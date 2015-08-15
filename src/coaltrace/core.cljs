(ns ^:figwheel-always coaltrace.core
    (:require [reagent.core :as r :refer [atom]]))

(enable-console-print!)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Helper Fns

(defn randhue []
  "make a hue from random color channel data"
  (str "#" (.toString (rand-int 16rFFFFFF) 16)))

(defn randnum []
  "get random num between -1 and 1"
  (let [ran (.random js/Math)
        plusminus (if (>  0.5 (.random js/Math)) 1 -1)]
    (* ran plusminus)))

(defn randnums-unique [max n]
  "get n unique ints less than or equal to max"
  (if (> n max)
    (into [] (range 0 max))
    (loop [numbers (take n (iterate rand-int (dec max)))]
      (if (= n (count (set numbers)))
        (into [] numbers)
        (recur (take n (iterate rand-int (dec max))))))))

(defn update-rand [coll f]
  "update a random member of a coll with function f"
  (let [index (rand-int (count coll))
        v (into [] coll)]
    (update-in v [index] f)))

(defn addvec [[x1 y1] [x2 y2]]
  "add two vectors"
  [(+ x1 x2) (+ y1 y2)])

(defn subtractvec [[x1 y1] [x2 y2]]
  "subtract two vectors"
  [(- x1 x2) (- y1 y2)])

(defn normalizevec [[x y]]
  "normalize a vector two vectors"
  (let [mag (.sqrt js/Math (+ (* x x) (* y y)))
        newx (/ x mag)
        newy (/ y mag)]
    [newx newy]))

(defn multiplyvec [[x y] s]
  "multiply vector by a scalar"
  [(* x s) (* y s)])

(defn constrain [val min max]
  "force the velocity vector betwen min and max vals"
  (cond (>= val max) max
        (<= val min) min
        :else val ))

(defn distance [ [x1 y1] [x2 y2] ]
  "distance of two points"
  (let [xdif (- x2 x1)
        ydif (- y2 y1)]
    (.sqrt js/Math (+ (* xdif xdif)) (* ydif ydif))))

(defn coulomb [d charge]
  "calculate coulomb force"
  (if (> d 0)
    (/ (.sqrt js/Math charge) (.sqrt js/Math d))
    1000))

(defn poissonSample [l]
  (let [t (.exp js/Math (* -1 l))]
    (loop [k 0 p 1]
      (if (< p t)
        (- k 1)
        (recur  (+ k 1) (* p (.random js/Math)))))))

(defn coal-interval [k N gen pushback]
  "coalescent intervals calculation, returns number of pixels in interval of k lineages
   every frame each member of trace is decremented by PUSHBACK pixels

   Moran model with overlapping generations, k concurrent lineages
  "
  (let [;Moran model with overlapping generations, k concurrent lineages
        cof (/ 2 (* k (- k 1)))
        generations (* cof (* 0.5 N))
        ;converting from generations to frames
        frames (* generations gen)
        ;converting from frames to pixels
        mod (* frames pushback)]
    mod))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Individuals and Their Operations

(defn reset-indv [i]
  (assoc i :vel [0 0] :acc [0 0]))

(defn replicate-indv [i]
  (let [[x y ] (:loc i)
        newx (+ x (randnum))
        newy (+ y (randnum))]
    (assoc i :loc [newx newy] :vel [0 0] :acc [0 0] :radius 0.01 :growing true :dying false)))

(defn die-indv [i]
  (assoc i :dying true :growing false))

(defn update-location [i two-D? height baseline maxvel]
  (let [v1 (addvec (:vel i) (:acc i))
        v2 (let [[x y] v1] [(constrain x (- maxvel) maxvel) (constrain y (- maxvel) maxvel)])
        loc1 (addvec v2  (:loc i))
        loc2 (let [[x y] loc1]
               (if two-D?
                 [x y]
                 [x (- height baseline)]))]
    (assoc i :loc loc2 :vel [0 0] :acc [0 0])))

(defn update-radius [i max]
  (let [r (:radius i)
        r1 (if (:growing i) (+ 0.9 r) r )
        r2 (if (>= r1 max) (- r1 0.4) r1)
        r3 (if (:dying i) (- r2 0.4) r2)]
    (assoc i :radius r3)))

(defn update-growing [i max]
  (if (>= (:radius i) max)
    (assoc i :growing false)
    i))

(defn update-trace [i maxtraces]
  "update the traces"
  (let [[x y] (:loc i)
        hue  (:hue i)
        t    (:trace i)
        newtrace (if (> (count t) maxtraces)
                   (conj (pop t) [x y hue])
                   (conj t [x y hue]))
        ]
    ;(if (> maxtraces (count t))
    ;  ;(assoc i :trace (conj newtrace (rest t)))
    ;  i
    ;  (assoc i :trace (conj newtrace t))
    (println t)
    (println (count t))
    (println maxtraces)
    (print newtrace)
      (assoc i :trace newtrace)))


(defn extend-trace [i]
  (let [[x y] (:loc i)
        hue  (:hue i)
        t    (:trace i)
        newtrace (-> (rest t)
                     (conj [x y hue]))]
    (assoc i :trace newtrace)))

(defn resettrace [i tracedepth]
  (let [[x y] (:loc i)
        hue (:hue i)
        newtrace (repeat tracedepth [x y hue])]
    (assoc i :trace newtrace)))

(defn mutate [i] (assoc i :hue (randhue)))

(defn update-all [i max two-D? height baseline maxvel]
  (-> i
      (update-growing max)
      (update-radius max)
      (update-location two-D? height baseline maxvel)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Functions that act on Population (vectors of Individuals)

(defn births
  "replicate n times"
  [pop n]
  ;{:pre [(> n 0)
  ;       (vector? pop)]
  ; :post [(vector? %)]}

  (let [p2 (shuffle pop)
        replicates (into [] (map replicate-indv (take n pop)))]
      (concat p2 replicates)))

(defn deaths
  "replicate n times"
  [pop n]
  ;{:pre [(> n 0)
  ;       (vector? pop)]
  ; :post [(vector? %)]}

  (let [p2 (shuffle pop)
        markdead (map die-indv (take n p2))]
    (concat markdead (drop n p2))))

(defn mutate-pop
  "mutate n times"
  [pop n]
  ;{:pre [(> n 0)
  ;       (vector? pop)]
  ; :post [(vector? %)] }

  (let [p2 (shuffle pop)
        mutated (map mutate (take n p2))]
    (concat mutated (drop n p2))))

(defn cleanup [population]
  "remove small radius"
  (filter #(<= 0 (:radius %))  population))

;; calc overlaps
(defn- overlap? [ind1 ind2]
  "check for overlap"
  (let [dist (distance (:loc ind1) (:loc ind2))
        r1 (:radius ind1)
        r2 (:radius ind2)
        overlap (- (+ r1 r2) dist)]
    (if (> overlap 0) true false)))

(defn exclude-points [population]
  "exclude overlaping populations"
  (let [p2 (into [] population)
        locs (for [i (range (count population))
                   j (range (count population))
                   :when (not= i j)] [i j])

        ;if overlapped, reset values
        overlapfn (fn [p loc]
                    (let [[i j] loc
                          p2 (into [] p)]
                      (if (overlap? (nth p2 i) (nth p2 j))
                        (-> p2 (update-in [i] reset-indv)
                            (update-in [j] reset-indv))
                        p2)))]

    (reduce overlapfn population locs)))

(defn exclude-walls [population & {:keys [width height]}]
  "keep points within wall boundaries"
  (let [exfn (fn [i] (let [[x y] (:loc i)
                           r (:radius i)
                           newx (constrain x (* r 2) (- width (* r 2)))
                           newy (constrain y (* r 2) (- height (* r 2)))]
                       (assoc i :loc [newx newy])))]
    (map exfn population)))

(defn repulsion-points [population charge]
  "update acceleration of an individual based on replusive forces from all others"
  (let [p2 (into [] population)
        accfn (fn [individual1 individual2]
                "obtain pairwise contribution to acceleration"
                (let [loc1 (:loc individual1)
                      loc2 (:loc individual2)
                      dist (distance loc1 loc2)]
                  (if (zero? dist)
                    nil
                    (let [force (coulomb dist charge)
                          normal (normalizevec (subtractvec loc2 loc1))
                          accvec (multiplyvec normal force)]
                      accvec))))

        update-ind (fn [i]
                     "add the pairwise contributions and update the acceleration"
                     (let [accvecs (map #(accfn % i) p2)
                           nonnil-accs (remove nil? accvecs)
                           newacc (reduce addvec (:acc i) nonnil-accs)]
                       (assoc i :acc newacc)))]
    (map update-ind p2)))

(defn repulsion-walls [population wallmultiplier width height charge]
  "add replusion from walls"
  (let [ wallrepfn (fn [ind]
                     "update points based on replusion from the wall."
                     (let [[x y] (:loc ind)
                           lwall [1 0]       rwall [-1 0]
                           twall [0 1]       bwall [0 -1]
                           ldist x           tdist y
                           rdist (- width x) bdist (- height y)
                           lacc (multiplyvec lwall (* wallmultiplier (coulomb ldist charge)))
                           racc (multiplyvec rwall (* wallmultiplier (coulomb rdist charge)))
                           tacc (multiplyvec twall (* wallmultiplier (coulomb tdist charge)))
                           bacc (multiplyvec bwall (* wallmultiplier (coulomb bdist charge)))
                           add-acc  (reduce addvec (:acc ind) [lacc racc tacc bacc])]
                       (assoc ind :acc add-acc)))]
    (map wallrepfn population)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Constants And App Control


(def constants
  (r/atom
    {:frame 1
     :two-d true
     :tracing true
     :statistics false
     :help  false
     :frate false
     :graybg true
     :branchcoloring true
     :looping true
     :dynamics true
     :mutation true
     :baseline 25
     :tracedepth 50
     :tracestep 18
     :pushback 0.75
     :gen 2
     :N 16
     :MU 0.1
     :charge 20
     :MAXVEL  2
     :MAXRAD  6
     :wallmultiplier 40
     :width  400
     :height 400}))

(def population-atom (r/atom []))

(defn make-individual [& {:keys [loc vel acc hue dying mutation trace growing radius]
                          :or   {loc [0 0]     vel [0 0]       acc [0 0]     hue  ""
                                 dying false   mutation false  trace #queue []      growing true
                                 radius 0.001 }}]
  {:loc  loc
   :vel vel
   :acc acc
   :hue hue
   :mutation mutation
   :trace trace
   :growing growing
   :dying dying
   :radius radius})

(defn newpopulation [n w h]
  "new population"
  (for [ _ (range n)]
    (let [x (rand-int w)
          y (rand-int h)]
      (make-individual :loc [x y] :hue "#22dd6f"))))

(let [data @constants
      popsize (:N data)
      width   (:width data)
      height  (:height data)
      population (into [] (newpopulation popsize width height)) ]
  (reset! population-atom population))

(defn update-population [pop]
  (let [dynamics       (:dynamics @constants)
        mutation       (:mutation @constants)
        gen            (:gen @constants)
        N              (:N @constants)
        MU             (:MU @constants)
        charge         (:charge @constants)
        wallmultiplier (:wallmultiplier @constants)
        width          (:width @constants)
        height         (:height @constants)
        two-D?         (:two-d @constants)
        rad-max        (:MAXRAD @constants)
        baseline       (:baseline 25)
        maxvel         (:MAXVEL @constants)
        births-deaths  (let [popBD (/ 1 (* gen N))] (poissonSample popBD))
        mutations      (let [popMU (/ 1 (* gen N MU))] (poissonSample popMU))
        p1 (if (and dynamics (> births-deaths 0))
             (-> pop (births births-deaths)
                     (deaths births-deaths))
             pop)
        p2 (if (and mutation (> mutations 0))
             (mutate-pop p1 mutations)
             p1)

        repulsed1 (repulsion-points p2 charge)
        repulsed2 (repulsion-walls repulsed1 wallmultiplier width height charge)
        updated (map #(update-all % rad-max two-D? height baseline maxvel) repulsed2)
        exclude1 (exclude-points updated)
        exclude2 (exclude-walls exclude1 :width width :height height)
        clean    (cleanup exclude2)
        updated-trace (map #(update-trace % (:tracedepth @constants)) clean)]
    updated-trace))


(defn circle-rep [individual]
  "circle representation"
  (let [[x y] (:loc individual)
        r (:radius individual)
        hue (:hue individual)]
    [:circle {:cx x :cy y :r r :fill hue :stroke "black"}] ))

(defn traces [individual]
  "trace representation"
  (let [traces (:trace individual)]
    (for [[t1 t2] (partition 2 traces)]
      (let [[x1 y1 h1] t1
            [x2 y2 h2] t2]
      [:line {:x1 x1  :y1 y1 :x2 x2 :y2 y2 :stroke h1 :stroke-width 5}]))))

(defn render-population []
  (fn []
    (js/setTimeout #(swap! population-atom update-population) 5)
    [:svg {:width  (:width  @constants)
           :height (:height @constants)}
     ;add traces
     (if (:tracing @constants)
       (for [ind @population-atom] (traces ind))
       nil)
     ;add circles
     (for [ind @population-atom]  (circle-rep ind))
     ]))

(r/render [render-population] (. js/document (getElementById "app")))
