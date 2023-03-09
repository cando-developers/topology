(in-package :topology)

;;; Copied from kinematics badgeom.lisp


(defun bad-carbon-geometry-p (aggregate)
  ;; Check that sp3 carbons have good geometry
  (let ((member-3-rings (let ((ht (make-hash-table)))
                          (loop for ring in chem:*current-rings*
                                when (= (length ring) 3)
                                  do (progn
                                       #+(or)(format t "ring3 ~a~%" ring)
                                       (loop for atm in ring
                                             do (setf (gethash atm ht) ring))))
                          ht))
        (member-4-rings (let ((ht (make-hash-table)))
                          (loop for ring in chem:*current-rings*
                                when (= (length ring) 4)
                                  do (progn
                                       (loop for atm in ring
                                             do (setf (gethash atm ht) ring))))
                          ht)))
    (labels
        ((calculate-angle-degrees (va vx vb)
           (/ (geom:calculate-angle va vx vb) 0.0174533))
         (bad-sp3-angle-p (angle n1 n2 n3)
           (when (or (< angle 90.0) (> angle 130.0))
             (return-from bad-carbon-geometry-p (format nil "bad sp3 angle(90.0.lt.angle.lt.130.0) ~7,2f between ~a-~a-~a"
                                                        angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3)))))
         (bad-sp3-angle-in-member-3-ring-p (angle in-member-3-ring n1 n2 n3)
           #+(or)(format t "bad-sp3-angle-in-member-3-ring-p in-member-3-ring ~a n1: ~a n2: ~a n3: ~a~%" in-member-3-ring n1 n2 n3)
           (if (and in-member-3-ring
                    (eq in-member-3-ring (gethash n1 member-3-rings))
                    (eq in-member-3-ring (gethash n3 member-3-rings)))
               (progn
                 (when (or (< angle 40.0) (> angle 90.0))
                   (return-from bad-carbon-geometry-p (format nil "member 3 ring bad sp3 angle(40.0.lt.angle.lt.90.0) ~7,2f between ~a-~a-~a"
                                                              angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3)))))
               (progn
                 (when (or (< angle 90.0) (> angle 130.0))
                   (return-from bad-carbon-geometry-p (format nil "member 3 ring bad sp3 angle(90.0.lt.angle.lt.130.0) ~7,2f between ~a-~a-~a"
                                                              angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3)))))))
         (bad-sp3-angle-in-member-4-ring-p (angle in-member-4-ring n1 n2 n3)
           (if (and in-member-4-ring
                    (eq in-member-4-ring (gethash n1 member-4-rings))
                    (eq in-member-4-ring (gethash n3 member-4-rings)))
               (when (or (< angle 60.0) (> angle 110.0))
                 (return-from bad-carbon-geometry-p (format nil "member 4 ring bad sp3 angle(60.0.lt.angle.lt.110.0) ~7,2f between ~a-~a-~a"
                                                            angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3))))
               (when (or (< angle 90.0) (> angle 130.0))
                 (return-from bad-carbon-geometry-p (format nil "member 4 ring bad sp3 angle(90.0.lt.angle.lt.130.0) ~7,2f between ~a-~a-~a"
                                                            angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3))))))
         (bad-sp2-angle-p (angle n1 n2 n3)
           (when (or (< angle 100.0) (> angle 140.0))
             (return-from bad-carbon-geometry-p (format nil "bad sp2 angle(100.0.lt.angle.lt.140.0) ~7,2f between ~a-~a-~a"
                                                        angle (chem:get-name n1) (chem:get-name n2) (chem:get-name n3)))))
         (bad-angles-p (aggregate)
           (chem:do-atoms (atm aggregate)
             (let ((vx (chem:get-position atm))
                   (in-member-3-ring (gethash atm member-3-rings))
                   (in-member-4-ring (gethash atm member-4-rings))
                   (bonds (chem:bonds-as-list atm)))
               #+(or)(format t "atm ~a in-member-3-ring ~a (length bonds) ~a~%" atm in-member-3-ring (length bonds))
               (cond
                 ((and (= (length bonds) 4) in-member-3-ring)
                  (let* ((aa (chem:bond/get-other-atom (first bonds) atm))
                         (ab (chem:bond/get-other-atom (second bonds) atm))
                         (ac (chem:bond/get-other-atom (third bonds) atm))
                         (ad (chem:bond/get-other-atom (fourth bonds) atm))
                         (va (chem:get-position aa))
                         (vb (chem:get-position ab))
                         (vc (chem:get-position ac))
                         (vd (chem:get-position ad))
                         (axb (calculate-angle-degrees va vx vb))
                         (axc (calculate-angle-degrees va vx vc))
                         (axd (calculate-angle-degrees va vx vd))
                         (bxc (calculate-angle-degrees vb vx vc))
                         (bxd (calculate-angle-degrees vb vx vd))
                         (cxd (calculate-angle-degrees vc vx vd)))
                    #+(or)(progn
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ab) axb)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ac) axc)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ad) axd)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ac) bxc)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ad) bxd)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ac) (chem:get-name atm) (chem:get-name ad) cxd))
                    (when (or (bad-sp3-angle-in-member-3-ring-p axb in-member-3-ring aa atm ab)
                              (bad-sp3-angle-in-member-3-ring-p axc in-member-3-ring aa atm ac)
                              (bad-sp3-angle-in-member-3-ring-p axd in-member-3-ring aa atm ad)
                              (bad-sp3-angle-in-member-3-ring-p bxc in-member-3-ring ab atm ac)
                              (bad-sp3-angle-in-member-3-ring-p bxd in-member-3-ring ab atm ad)
                              (bad-sp3-angle-in-member-3-ring-p cxd in-member-3-ring ac atm ad))
                      (return-from bad-angles-p t))
                    ))
                 ((and (= (length bonds) 4) in-member-4-ring)
                  (let* ((aa (chem:bond/get-other-atom (first bonds) atm))
                         (ab (chem:bond/get-other-atom (second bonds) atm))
                         (ac (chem:bond/get-other-atom (third bonds) atm))
                         (ad (chem:bond/get-other-atom (fourth bonds) atm))
                         (va (chem:get-position aa))
                         (vb (chem:get-position ab))
                         (vc (chem:get-position ac))
                         (vd (chem:get-position ad))
                         (axb (calculate-angle-degrees va vx vb))
                         (axc (calculate-angle-degrees va vx vc))
                         (axd (calculate-angle-degrees va vx vd))
                         (bxc (calculate-angle-degrees vb vx vc))
                         (bxd (calculate-angle-degrees vb vx vd))
                         (cxd (calculate-angle-degrees vc vx vd)))
                    #+(or)(progn
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ab) axb)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ac) axc)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ad) axd)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ac) bxc)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ad) bxd)
                            (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ac) (chem:get-name atm) (chem:get-name ad) cxd))
                    (when (or (bad-sp3-angle-in-member-4-ring-p axb in-member-4-ring aa atm ab)
                              (bad-sp3-angle-in-member-4-ring-p axc in-member-4-ring aa atm ac)
                              (bad-sp3-angle-in-member-4-ring-p axd in-member-4-ring aa atm ad)
                              (bad-sp3-angle-in-member-4-ring-p bxc in-member-4-ring ab atm ac)
                              (bad-sp3-angle-in-member-4-ring-p bxd in-member-4-ring ab atm ad)
                              (bad-sp3-angle-in-member-4-ring-p cxd in-member-4-ring ac atm ad))
                      (return-from bad-angles-p t))))
                 ((= (length bonds) 4)
                  (let* ((aa (chem:bond/get-other-atom (first bonds) atm))
                         (ab (chem:bond/get-other-atom (second bonds) atm))
                         (ac (chem:bond/get-other-atom (third bonds) atm))
                         (ad (chem:bond/get-other-atom (fourth bonds) atm))
                         (va (chem:get-position aa))
                         (vb (chem:get-position ab))
                         (vc (chem:get-position ac))
                         (vd (chem:get-position ad))
                         (axb (calculate-angle-degrees va vx vb))
                         (axc (calculate-angle-degrees va vx vc))
                         (axd (calculate-angle-degrees va vx vd))
                         (bxc (calculate-angle-degrees vb vx vc))
                         (bxd (calculate-angle-degrees vb vx vd))
                         (cxd (calculate-angle-degrees vc vx vd)))
                    (or (bad-sp3-angle-p axb aa atm ab)
                        (bad-sp3-angle-p axc aa atm ac)
                        (bad-sp3-angle-p axd aa atm ad)
                        (bad-sp3-angle-p bxc ab atm ac)
                        (bad-sp3-angle-p bxd ab atm ad)
                        (bad-sp3-angle-p cxd ac atm ad))))
                 ((= (length bonds) 3)
                  (let* ((aa (chem:bond/get-other-atom (first bonds) atm))
                         (ab (chem:bond/get-other-atom (second bonds) atm))
                         (ac (chem:bond/get-other-atom (third bonds) atm))
                         (va (chem:get-position aa))
                         (vb (chem:get-position ab))
                         (vc (chem:get-position ac))
                         (axb (calculate-angle-degrees va vx vb))
                         (axc (calculate-angle-degrees va vx vc))
                         (bxc (calculate-angle-degrees vb vx vc)))
                    (or (bad-sp2-angle-p axb aa atm ab)
                        (bad-sp2-angle-p axc aa atm ac)
                        (bad-sp2-angle-p bxc ab atm ac)))))))))
      (or (bad-angles-p aggregate)
          ))))

(defparameter *dkp* (chem:make-chem-info-graph (chem:compile-smarts "[C&H1:0]1-[N:1](-[C:2])-[C:3](=[O:4])-[C&H0:5]-[N:6](-[C&H2:7])-[C:8]1(=[O:9])")))

(defun bad-dkp-geometry-p (agg)
  (let* ((mol (chem:content-at agg 0))
         (mol-graph (chem:make-molecule-graph-from-molecule mol))
         (matches (chem:boost-graph-vf2 *dkp* mol-graph)))
    (loop for match in matches
          for c1 = (elt match 2)
          for n1 = (elt match 1)
          for co1 = (elt match 3)
          for oc1 = (elt match 4)
          for c1n = (chem:get-name c1)
          for n1n = (chem:get-name n1)
          for co1n = (chem:get-name co1)
          for oc1n = (chem:get-name oc1)
          for dihed1 = (abs (/ (chem:calculate-dihedral-for-atoms c1 n1 co1 oc1) 0.0174533))
          for c2 = (elt match 7)
          for n2 = (elt match 6)
          for co2 = (elt match 8)
          for oc2 = (elt match 9)
          for c2n = (chem:get-name c2)
          for n2n = (chem:get-name n2)
          for co2n = (chem:get-name co2)
          for oc2n = (chem:get-name oc2)
          for dihed2 = (abs (/ (chem:calculate-dihedral-for-atoms c2 n2 co2 oc2) 0.0174533))
          when (or (> dihed1 45.0)
                   (> dihed2 45.0))
            do (return-from bad-dkp-geometry-p
                 (format nil "bad-aromatic-geometry-p ~a-~a-~a-~a ~a>45.0  or ~a-~a-~a-~a ~a>45.0"
                         c1n n1n co1n oc1n (/ dihed1 0.0174533)
                         c2n n2n co2n oc2n (/ dihed2 0.0174533)
                         ))))
  nil)

(defparameter *planar-aromatic* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1:[C&H1:1](~[#1:2]):[C&H1:3](-[#1:4]):[*:5]:[*]:[*]1")))

(defun bad-aromatic-geometry-p (agg)
  (unless chem:*current-rings* (error "The chem:*current-rings* dynamic variable must be defined - use (chem:identify-rings matter)"))
  (chem:with-aromaticity-information (agg :am1bcc)
    (let* ((mol (chem:content-at agg 0))
           (mol-graph (chem:make-molecule-graph-from-molecule mol))
           (matches (chem:boost-graph-vf2 *planar-aromatic* mol-graph)))
      (loop for match in matches
            for h1 = (elt match 2)
            for c1 = (elt match 1)
            for c2 = (elt match 3)
            for h2 = (elt match 4)
            for h1n = (chem:get-name h1)
            for c1n = (chem:get-name c1)
            for h2n = (chem:get-name h2)
            for c2n = (chem:get-name c2)
            for ph1 = (chem:get-position h1)
            for pc1 = (chem:get-position c1)
            for pc2 = (chem:get-position c2)
            for ph2 = (chem:get-position h2)
            for torsion = (geom:calculate-dihedral ph1 pc1 pc2 ph2)
            when (> (abs (/ torsion 0.0174533)) 30.0)
              do (return (format nil "bad-aromatic-geometry-p ~a-~a-~a-~a ~a" h1n c1n c2n h2n (/ torsion 0.0174533)))
            finally (return nil)))))


(defun bad-geometry-p (agg)
  "Test for different kinds of common bad geometry 
- return true if there was any bad geometry and nil if there wasn't"
  (or (let ((chem:*current-rings* (chem:identify-rings agg)))
        (when chem:*current-rings*
          (or (bad-carbon-geometry-p agg)
              (bad-aromatic-geometry-p agg)
              (bad-dkp-geometry-p agg))))))
