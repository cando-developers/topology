(in-package :foldamer)

(defun build-monomer-context-index-map (confs)
  (let ((monomer-context-index (make-hash-table :test 'equal))
        (monomer-contexts-vector (make-array (hash-table-count (topology:monomer-context-to-fragment-conformations confs)))))
    (with-hash-table-iterator (my-iterator (topology:monomer-context-to-fragment-conformations confs))
      (loop
        named indexes
        for index from 0
        do (multiple-value-bind (entry-p key value)
               (my-iterator)
             (declare (ignore value))
             (if entry-p
                 (progn
                   (setf (gethash key monomer-context-index) index)
                   (setf (aref monomer-contexts-vector index) key))
                 (return-from indexes)))))
    (values monomer-context-index monomer-contexts-vector)))

(defun angle-difference (b1 b2)
  (let ((diff (mod (- b2 b1) 360)))
    (if (< diff -180)
	(incf diff 360)
	(if (> diff 180)
	    (decf diff 360)
	    diff))))

(defun fragments-match-p (before-fragment after-fragment after-fragment-focus-monomer-name)
  "Return T for the time being"
  (let ((before-fragment-out-of-focus (gethash after-fragment-focus-monomer-name (topology:out-of-focus-internals before-fragment))))
    (unless before-fragment-out-of-focus
      (error "Could not find ~s in out-of-focus-internals of ~s"
             after-fragment-focus-monomer-name before-fragment))
    (loop for index below (length before-fragment-out-of-focus)
          for before-bonded-internal = (elt before-fragment-out-of-focus index)
          for after-bonded-internal = (elt (topology:internals after-fragment) index)
          do (progn
               (unless (eq (topology:name before-bonded-internal)
                           (topology:name after-bonded-internal))
                 (error "Atom names don't match ~s ~s" before-bonded-internal after-bonded-internal))
               (let* ((diha (/ (topology:dihedral before-bonded-internal) 0.0174533))
                      (dihb (/ (topology:dihedral after-bonded-internal) 0.0174533))
                      (delta-angle (abs (angle-difference diha dihb))))
                 (when (> delta-angle 20.0)
                   (return-from fragments-match-p nil)))))
    t))

(defun match-conformations (before-monomer-context before-fragment-conformations
                            after-monomer-context after-fragment-conformations
                            fragment-conformations-map verbose)
  (loop with before-monomer-context-index = (gethash before-monomer-context (topology:monomer-context-index-map fragment-conformations-map))
        with after-monomer-context-index = (gethash after-monomer-context (topology:monomer-context-index-map fragment-conformations-map))
        with before-fragments = (topology:fragments before-fragment-conformations)
        with after-fragments =  (topology:fragments after-fragment-conformations)
        with before-fragment-match-table = (make-array (length before-fragments))
        with after-fragment-focus-monomer-name = (topology:focus-monomer-name after-fragment-conformations)
        for before-fragment across before-fragments
        for before-fragment-index from 0
        for after-fragment-match-vector = (loop named build-or-reuse-match
                                                with new-after-fragment-match-vector = (make-array 16 :element-type 'fixnum :adjustable t :fill-pointer 0)
                                                for after-fragment across after-fragments
                                                for after-fragment-index from 0
                                                for match = (fragments-match-p before-fragment after-fragment after-fragment-focus-monomer-name)
                                                when match
                                                  do (vector-push-extend after-fragment-index new-after-fragment-match-vector)
                                                finally (loop named reuse-results
                                                              for ii from 0 below before-fragment-index
                                                              for previous-after-fragment-match-vector = (elt before-fragment-match-table ii)
                                                              when (equalp previous-after-fragment-match-vector new-after-fragment-match-vector)
                                                                do (return-from build-or-reuse-match previous-after-fragment-match-vector)
                                                              finally (return-from build-or-reuse-match new-after-fragment-match-vector)))
        do (when (and verbose (= (length after-fragment-match-vector) 0))
             (format t "Empty after-fragment-match-vector for ~a -> ~a : ~a~%" before-monomer-context after-monomer-context before-fragment-index)
             (let* ((fragment (elt (topology:fragments before-fragment-conformations) before-fragment-index))
                    (out-of-focus-internals (gethash after-fragment-focus-monomer-name (topology:out-of-focus-internals fragment))))
               (format t "~14a~{(~{~3a ~7,1f~})~}~%"
                       "Before"
                       (loop for internal across out-of-focus-internals
                             collect (list (topology:name internal) (/ (topology:dihedral internal) 0.0174533))))
               (loop for index below (length (topology:fragments after-fragment-conformations))
                     for after-fragment = (elt (topology:fragments after-fragment-conformations) index)
                     do (format t "~10a~4a~{(~{~3a ~7,1f~})~}~%"
                                "After" (format nil "~-4d" index)
                                (loop for ii below 5
                                      for internal = (elt (topology:internals after-fragment) ii)
                                      collect (list (topology:name internal) (/ (topology:dihedral internal) 0.0174533)))))
               #+(or)(break "out-of-focus-internals ~a" out-of-focus-internals)
               ))
        if (= (length after-fragment-match-vector) 0)
          do (pushnew before-fragment-index (gethash (topology:make-fragment-match-key
                                                      :before-monomer-context-index before-monomer-context-index
                                                      :after-monomer-context-index after-monomer-context-index)
                                                     (topology:missing-fragment-matches fragment-conformations-map)
                                                     nil))
        else
          do (let ((fragment-match-key (topology:make-fragment-match-key
                                        :before-monomer-context-index before-monomer-context-index
                                        :after-monomer-context-index after-monomer-context-index)))
               (let ((total-after-fragment-match-vector (gethash fragment-match-key (topology:fragment-matches fragment-conformations-map))))
                 (cond
                   ((null total-after-fragment-match-vector)
                    (setf total-after-fragment-match-vector (make-array (1+ before-fragment-index) :adjustable t)))
                   ((<= (length total-after-fragment-match-vector) before-fragment-index)
                    (adjust-array total-after-fragment-match-vector (1+ before-fragment-index))))
                 (setf (aref total-after-fragment-match-vector before-fragment-index)
                       after-fragment-match-vector)
                 (setf (gethash fragment-match-key (topology:fragment-matches fragment-conformations-map)) total-after-fragment-match-vector))))
  fragment-conformations-map)


(defun matching-oligomers-p (before-monomer-context before-oligomer before-focus-monomer before-focus-monomer-name
                             after-monomer-context after-oligomer after-focus-monomer after-focus-monomer-name)
  (declare (ignore before-monomer-context after-monomer-context))
  (let ((after-iterator (topology:directional-coupling-iterator after-oligomer))
        after-coupling after-source-monomer-name after-target-monomer-name
        before-coupling before-source-monomer-name before-target-monomer-name)
    ;; Find the in-coupling for the after-oligomer from the after-focus-monomer
    (loop named after-loop
          do (multiple-value-setq (after-coupling after-source-monomer-name after-target-monomer-name) (funcall after-iterator))
          do (unless after-coupling (return-from matching-oligomers-p nil))
          do (when (eq (topology:target-monomer after-coupling) after-focus-monomer)
               (return-from after-loop nil)))
    (let ((before-iterator (topology:directional-coupling-iterator before-oligomer)))
      (loop named before-loop
            do (multiple-value-setq (before-coupling before-source-monomer-name before-target-monomer-name) (funcall before-iterator))
            do (unless before-coupling (return-from matching-oligomers-p nil))
            do (when (and (eq (topology:source-monomer before-coupling) before-focus-monomer)
                          (eq (topology:source-plug-name before-coupling) (topology:source-plug-name after-coupling)))
                 (return-from before-loop nil)))
      (let ((before-target-monomer-name (topology:oligomer-monomer-name-for-monomer before-oligomer (topology:target-monomer before-coupling)))
            (after-source-monomer-name (topology:oligomer-monomer-name-for-monomer after-oligomer (topology:source-monomer after-coupling))))
        (and (eq before-target-monomer-name after-focus-monomer-name)
             (eq after-source-monomer-name before-focus-monomer-name))))))


(defclass monomer-context-pair ()
  ((before-monomer-context :initarg :before-monomer-context :accessor before-monomer-context)
   (after-monomer-context :initarg :after-monomer-context :accessor after-monomer-context)))

(defun over-matching-oligomers (before-training-oligomer-space
                                after-training-oligomer-space)
  (let ((before-iterator (oligomer-monomer-context-focus-monomer-iterator before-training-oligomer-space))
        result)
    (loop named before-loop
          do (multiple-value-bind (number-remaining before-oligomer before-monomer-context before-focus-monomer before-focus-monomer-name)
                 (funcall before-iterator)
               (unless number-remaining (return-from before-loop nil))
               (when (and number-remaining (> number-remaining 0))
                 (let ((after-iterator (oligomer-monomer-context-focus-monomer-iterator after-training-oligomer-space)))
                   (loop named after-loop
                         do (multiple-value-bind (number-remaining after-oligomer after-monomer-context after-focus-monomer after-focus-monomer-name)
                                (funcall after-iterator)
                              (unless number-remaining (return-from after-loop nil))
                              (when (and number-remaining (> number-remaining 0))
                                (when (matching-oligomers-p before-monomer-context before-oligomer before-focus-monomer before-focus-monomer-name
                                                            after-monomer-context after-oligomer after-focus-monomer after-focus-monomer-name)
                                  (push (make-instance 'monomer-context-pair
                                                       :before-monomer-context before-monomer-context
                                                       :after-monomer-context after-monomer-context)
                                        result)
                                  #+(or)(format t "Do something with: ~a -> ~a~%" before-monomer-context after-monomer-context)))))))))
    result))

(defun all-matching-monomer-contexts (foldamer)
  (loop named indexes
        for index from 0
        for training-oligomer-space in (training-oligomer-spaces foldamer)
        for oligomer-space = (oligomer-space training-oligomer-space)
        for focus-monomer = (focus-monomer training-oligomer-space)
        for focus-monomer-id = (topology:id focus-monomer)
        for in-coupling = (loop named find-in-coupling
                                for coupling across (topology:couplings oligomer-space)
                                when (eq (topology:target-monomer coupling) focus-monomer)
                                  do (return-from find-in-coupling coupling))
        when in-coupling
          append (let* ((in-monomer (topology:source-monomer in-coupling))
                        (in-monomer-id (topology:id in-monomer)))
                   (declare (ignore in-monomer-id))
                   (loop for in-training-oligomer-space in (training-oligomer-spaces foldamer)
                         for in-oligomer-space = (oligomer-space in-training-oligomer-space)
                         for in-focus-monomer = (focus-monomer in-training-oligomer-space)
                         for in-focus-monomer-id = (topology:id in-focus-monomer)
                         append (loop named find-out-coupling
                                      for coupling across (topology:couplings in-oligomer-space)
                                      for out-monomer = (topology:target-monomer coupling)
                                      for out-monomer-id = (topology:id out-monomer)
                                      when (and (eq (topology:id (topology:source-monomer coupling)) in-focus-monomer-id)
                                                (eq (topology:source-plug-name coupling) (topology:source-plug-name in-coupling))
                                                (eq (topology:target-plug-name coupling) (topology:target-plug-name in-coupling))
                                                (eq (topology:id (topology:source-monomer coupling))
                                                    (topology:id (topology:source-monomer in-coupling)))
                                                (eq (topology:id (topology:target-monomer coupling))
                                                    (topology:id (topology:target-monomer in-coupling))))
                                        append (over-matching-oligomers in-training-oligomer-space
                                                                        training-oligomer-space)
                                      #+(or)(all-matching-momomer-contexts-in-matching-training-oligomer-space-pairs
                                             in-training-oligomer-space
                                             training-oligomer-space))))))

(defun optimize-fragment-conformations-map (simple-fragment-conformations-map foldamer verbose)
  (let ((all-matching-monomer-contexts (all-matching-monomer-contexts foldamer)))
    (multiple-value-bind (monomer-context-index-map monomer-contexts-vector)
        (build-monomer-context-index-map simple-fragment-conformations-map)
      (let ((fragment-conformations-map
              (make-instance 'topology:matched-fragment-conformations-map
                             :monomer-context-to-fragment-conformations
                             (topology:monomer-context-to-fragment-conformations simple-fragment-conformations-map)
                             :monomer-context-index-map monomer-context-index-map
                             :monomer-contexts-vector monomer-contexts-vector))
            (num-pairs 0))
        (loop for monomer-context-pair in all-matching-monomer-contexts
              for before-monomer-context = (before-monomer-context monomer-context-pair)
              for after-monomer-context = (after-monomer-context monomer-context-pair)
              with monomer-context-to-fragment-conformations = (topology:monomer-context-to-fragment-conformations fragment-conformations-map)
              for before-monomer-context-index = (gethash before-monomer-context monomer-context-index-map)
              for before-fragment-conformations = (gethash before-monomer-context monomer-context-to-fragment-conformations) 
              for after-monomer-context-index = (gethash after-monomer-context monomer-context-index-map)
              for after-fragment-conformations = (gethash after-monomer-context monomer-context-to-fragment-conformations)
              do (match-conformations before-monomer-context before-fragment-conformations
                                      after-monomer-context after-fragment-conformations
                                      fragment-conformations-map verbose)
              do (incf num-pairs))
        (values fragment-conformations-map num-pairs)))))





(defun analyze-missing-conformations (confs)
  (let ((*print-pretty* nil))
    (maphash (lambda (key value)
               (let* ((before-monomer-context-index (topology:missing-fragment-match-key-before-monomer-context-index key))
                      (after-monomer-context-index (topology:missing-fragment-match-key-after-monomer-context-index key))
                      (before-monomer-context (elt (topology:monomer-contexts-vector confs) before-monomer-context-index))
                      (after-monomer-context (elt (topology:monomer-contexts-vector confs) after-monomer-context-index))
                      (before-fragment (gethash before-monomer-context (topology:monomer-context-to-fragment-conformations confs)))
                      (after-fragment (gethash after-monomer-context (topology:monomer-context-to-fragment-conformations confs))))
                 (format t "~a -> ~a : ~s~%" before-monomer-context after-monomer-context value)
                 (loop for internals-index in value
                       for fragment = (elt (topology:fragments before-fragment) internals-index)
                       do (format t "~a~%" (topology:out-of-focus-internals fragment))
                       do (break "examine ~a ~a" fragment after-fragment)
                       )
                 ))
             (topology:missing-fragment-matches confs))))

