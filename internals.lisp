(in-package :topology)

(defclass internal ()
  ((name :initarg :name :accessor name)
   ))

(cando:make-class-save-load
 internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass jump-internal (internal)
  ())

(cando:make-class-save-load jump-internal)

(defclass bonded-internal (internal)
  ((bond :initarg :bond :accessor bond)
   (angle :initarg :angle :accessor angle)
   (dihedral :initarg :dihedral :accessor dihedral)))

(cando:make-class-save-load
 bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass complex-bonded-internal (bonded-internal)
  ())

(cando:make-class-save-load complex-bonded-internal)

(defclass out-of-focus-internal ()
  ((name :initarg :name :accessor name)
   (atres-name :initarg :atres-name :accessor atres-name)
   (p-name :initarg :p-name :accessor p-name)
   (p-atres-name :initarg :p-atres-name :accessor p-atres-name)
   (gp-name :initarg :gp-name :accessor gp-name)
   (gp-atres-name :initarg :gp-atres-name :accessor gp-atres-name)
   (ggp-name :initarg :ggp-name :accessor ggp-name)
   (ggp-atres-name :initarg :ggp-atres-name :accessor ggp-atres-name)
   (dihedral-rad :initarg :dihedral-rad :accessor dihedral-rad)))

(cando:make-class-save-load out-of-focus-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a(~a) ~a(~a) ~a(~a) ~a(~a) ~5,2f"
             (name obj)
             (atres-name obj)
             (p-name obj)
             (p-atres-name obj)
             (gp-name obj)
             (gp-atres-name obj)
             (ggp-name obj)
             (ggp-atres-name obj)
             (dihedral-rad obj)))))

(defclass fragment-internals ()
  ((index :initarg :index :accessor index)
   (internals :initarg :internals :accessor internals)
   (out-of-focus-internals :initarg :out-of-focus-internals :accessor out-of-focus-internals)))

(cando:make-class-save-load fragment-internals)

(defclass fragment-conformations ()
  ((focus-monomer-name :initarg :focus-monomer-name :accessor focus-monomer-name)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (total-count :initform 0 :initarg :total-count :accessor total-count)
   (fragments :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :fragments :accessor fragments)))

(cando:make-class-save-load fragment-conformations)

(defclass fragment-conformations-map ()
  ((monomer-context-to-fragment-conformations :initform (make-hash-table :test 'equal)
                                              :initarg :monomer-context-to-fragment-conformations
                                              :accessor monomer-context-to-fragment-conformations)))

(cando:make-class-save-load fragment-conformations-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))

(defclass matched-fragment-conformations-map (fragment-conformations-map)
  ((fragment-matches :initform (make-hash-table :test 'equalp)
                     :initarg :fragment-matches
                     :accessor fragment-matches)))

(cando:make-class-save-load matched-fragment-conformations-map)

(defun matched-fragment-conformations-summary (matched-fragment-conformations-map)
  (let ((total-fragment-conformations 0)
        (matching-fragment-conformations 0)
        (missing-fragment-conformations 0))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf total-fragment-conformations (length (fragments value))))
             (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf matching-fragment-conformations (length value)))
             (fragment-matches matched-fragment-conformations-map))
    (maphash (lambda (key value)
               (declare (ignorable key))
               (loop for val across value
                     when (= (length val) 0)
                       do (incf missing-fragment-conformations)))
             (fragment-matches matched-fragment-conformations-map))
    (let ((missing-monomer-contexts nil))
      (maphash (lambda (monomer-context fragment-conformations)
                 (block inner-search
                   (maphash (lambda (monomer-context-pair allowed-fragment-indices)
                              (when (or (string= (car monomer-context-pair) monomer-context)
                                        (string= (cdr monomer-context-pair) monomer-context))
                                (return-from inner-search nil)))
                            (fragment-matches matched-fragment-conformations-map))
                   (push monomer-context missing-monomer-contexts)))
               (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
      (values total-fragment-conformations matching-fragment-conformations missing-fragment-conformations missing-monomer-contexts))))

(defconstant +dihedral-threshold+ (* 10.0 0.0174533))

(defun similar-internals-p (frag1 frag2 &optional )
  (loop for frag1-int across (internals frag1)
        for frag2-int across (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((aa (- (dihedral frag1-int) (dihedral frag2-int)))
                    (aamod (- (mod (+ aa 180) 360) 180)))
               (when (> (abs aamod) +dihedral-threshold+)
                 (return-from similar-internals-p nil)))))
  t)

(defun seen-fragment-internals (fragment-conformations fragment-internals)
  (loop for seen-frag across (fragments fragment-conformations)
        when (similar-internals-p seen-frag fragment-internals)
          do (return-from seen-fragment-internals (index seen-frag)))
  nil)


(defun bad-fragment-internals (fragment-internals)
  (loop for internal across (internals fragment-internals)
        with previous-internal
        do (cond
             ((typep internal 'jump-internal) nil)
             ((> (bond internal) 3.0)
              (return-from bad-fragment-internals
                (if previous-internal
                    (format nil "For atom ~a to ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) (name previous-internal) (bond internal))
                    (format nil "For atom ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) (bond internal))))))
        do (setf previous-internal internal))
  nil)

(defun save-fragment-conformations (fragment-conformations filename)
  (cando:save-cando fragment-conformations filename))

(defun load-fragment-conformations (filename)
  (cando:load-cando filename))

(defun dump-fragment-internals (fragment-internals finternals)
  (format finternals "begin-conformation ~a~%" (index fragment-internals))
  (flet ((to-deg (rad)
           (/ rad 0.0174533)))
    (loop for internal across (topology:internals fragment-internals)
          do (cond
               ((typep internal 'topology:jump-internal)
                (format finternals "jump-joint ~a~%" (topology:name internal)))
               ((typep internal 'topology:complex-bonded-internal)
                (format finternals "complex-bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                        (topology:name internal)
                        (topology:bond internal)
                        (to-deg (topology:angle internal))
                        (to-deg (topology:dihedral internal))))
               ((typep internal 'topology:bonded-internal)
                (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                        (topology:name internal)
                        (topology:bond internal)
                        (to-deg (topology:angle internal))
                        (to-deg (topology:dihedral internal))))
               ))
    (format finternals "end-conformation~%")
    (maphash (lambda (key internals)
               (format finternals "out-of-focus following-monomer ~a~%" key)
               (loop for internal across internals
                     do (cond
                          ((typep internal 'topology:jump-internal)
                           (format finternals "jump-joint ~a~%" (topology:name internal)))
                          ((typep internal 'topology:complex-bonded-internal)
                           (format finternals "complex-bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                                   (topology:name internal)
                                   (topology:bond internal)
                                   (to-deg (topology:angle internal))
                                   (to-deg (topology:dihedral internal))))
                          ((typep internal 'topology:bonded-internal)
                           (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                                   (topology:name internal)
                                   (topology:bond internal)
                                   (to-deg (topology:angle internal))
                                   (to-deg (topology:dihedral internal))))
                          )))
             (topology:out-of-focus-internals fragment-internals))))


(defgeneric fill-joint-internals (joint internal))

(defmethod fill-joint-internals ((joint kin:jump-joint) (internal jump-internal))
  )


(defmethod fill-joint-internals ((joint kin:bonded-joint) (internal bonded-internal))
  (kin:set-distance joint (bond internal))
  (kin:set-theta joint (angle internal))
  (kin:set-phi joint (dihedral internal))
  )
