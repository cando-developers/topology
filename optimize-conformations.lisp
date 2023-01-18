(in-package :foldamer)

(defun build-monomer-context-index-map (confs)
  (let ((monomer-context-index (make-hash-table :test 'equal)))
    (with-hash-table-iterator (my-iterator (topology:monomer-context-to-fragment-conformations confs))
      (loop
        named indexes
        for index from 0
        do (multiple-value-bind (entry-p key value)
               (my-iterator)
             (declare (ignore value))
             (if entry-p
                 (setf (gethash key monomer-context-index) index)
                 (return-from indexes)))))
    monomer-context-index))


(defclass training-oligomer-space-pair ()
  ((before-training-oligomer-space :initarg :before-training-oligomer-space :accessor before-training-oligomer-space)
   (after-training-oligomer-space :initarg :after-training-oligomer-space :accessor after-training-oligomer-space)))


(defun all-training-oligomer-space-pairs (foldamer)
  (let ((pairs nil))
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
            do (let* ((in-monomer (topology:source-monomer in-coupling))
                      (in-monomer-id (topology:id in-monomer)))
                 (declare (ignore in-monomer-id))
                 (loop for in-training-oligomer-space in (training-oligomer-spaces foldamer)
                       for in-oligomer-space = (oligomer-space in-training-oligomer-space)
                       for in-focus-monomer = (focus-monomer in-training-oligomer-space)
                       for in-focus-monomer-id = (topology:id in-focus-monomer)
                       for out-coupling = (loop named find-out-coupling
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
                                                  do (let* ((pair (make-instance 'training-oligomer-space-pair
                                                                                 :before-training-oligomer-space in-training-oligomer-space
                                                                                 :after-training-oligomer-space training-oligomer-space)))
                                                       (push pair pairs))))))
    pairs))


(defun fragments-match-p (before-fragment after-fragment)
  "Return T for the time being"
  (break "Check before-fragment ~a after-fragment ~a" before-fragment after-fragment))


(defun match-conformations (before-monomer-context before-fragment-conformations
                            after-monomer-context after-fragment-conformations
                            fragment-conformations-map)
  (loop with before-monomer-context-index = (gethash before-monomer-context (topology:monomer-context-index-map fragment-conformations-map))
        with after-monomer-context-index = (gethash after-monomer-context (topology:monomer-context-index-map fragment-conformations-map))
        with before-fragments = (topology:fragments before-fragment-conformations)
        with after-fragments =  (topology:fragments after-fragment-conformations)
        with before-fragment-match-table = (make-array (length before-fragments))
        for before-fragment in before-fragments
        for before-fragment-index from 0
        for after-fragment-match-vector = (loop named build-or-reuse-match
                                                with new-after-fragment-match-vector = (make-array 16 :element-type 'fixnum :adjustable t :fill-pointer 0)
                                                for after-fragment in after-fragments
                                                for after-fragment-index from 0
                                                for match = (fragments-match-p before-fragment after-fragment)
                                                when match
                                                  do (vector-push-extend after-fragment-index new-after-fragment-match-vector)
                                                finally (loop named reuse-results
                                                              for ii from 0 below before-fragment-index
                                                              for previous-after-fragment-match-vector = (elt before-fragment-match-table ii)
                                                              when (equalp previous-after-fragment-match-vector new-after-fragment-match-vector)
                                                                do (return-from build-or-reuse-match previous-after-fragment-match-vector)
                                                              finally (return-from build-or-reuse-match new-after-fragment-match-vector)))
        if (= (length after-fragment-match-vector) 0)
          do (push (make-instance 'topology:missing-fragment-match
                                  :before-monomer-context before-monomer-context
                                  :after-monomer-context after-monomer-context
                                  :before-fragment-conformation-index before-fragment-index)
                   (topology:missing-fragment-matches fragment-conformations-map))
        else
          do (let ((fragment-match-key (topology:make-fragment-match-key
                                        :before-monomer-context-index before-monomer-context-index
                                        :after-monomer-context-index after-monomer-context-index
                                        :before-fragment-conformation-index before-fragment-index)))
               (setf (gethash fragment-match-key (topology:fragment-matches fragment-conformations-map))
                     after-fragment-match-vector)))
  fragment-conformations-map)

(defun optimize-fragment-conformations-map (simple-fragment-conformations-map foldamer)
  (let* ((num-pairs 0)
         (training-oligomer-space-to-monomer-conformations (make-hash-table))
         (monomer-context-index-map (build-monomer-context-index-map simple-fragment-conformations-map))
         (fragment-conformations-map
           (make-instance 'topology:matched-fragment-conformations-map
                          :monomer-context-to-fragment-conformations
                          (topology:monomer-context-to-fragment-conformations simple-fragment-conformations-map)
                          :monomer-context-index-map monomer-context-index-map)))
    (loop for training-oligomer-space in (training-oligomer-spaces foldamer)
          for monomer-contexts = (foldamer:all-monomer-contexts-in-training-oligomer-space training-oligomer-space)
          do (setf (gethash training-oligomer-space training-oligomer-space-to-monomer-conformations) monomer-contexts))
    (loop for pair in (all-training-oligomer-space-pairs foldamer)
          for before-monomer-contexts = (all-monomer-contexts-in-training-oligomer-space (before-training-oligomer-space pair))
          for after-monomer-contexts = (all-monomer-contexts-in-training-oligomer-space (after-training-oligomer-space pair))
          with monomer-context-to-fragment-conformations = (topology:monomer-context-to-fragment-conformations fragment-conformations-map)
          do (loop for before-monomer-context in before-monomer-contexts
                   for before-monomer-context-index = (gethash before-monomer-context monomer-context-index-map)
                   for before-fragment-conformations = (gethash before-monomer-context monomer-context-to-fragment-conformations) 
                   do (loop for after-monomer-context in after-monomer-contexts
                            for after-monomer-context-index = (gethash after-monomer-context monomer-context-index-map)
                            for after-fragment-conformations = (gethash after-monomer-context monomer-context-to-fragment-conformations)
                            do (match-conformations before-monomer-context before-fragment-conformations
                                                    after-monomer-context after-fragment-conformations
                                                    fragment-conformations-map)
                            do (incf num-pairs)
                            )
                   )
          )
    (values fragment-conformations-map num-pairs)
    ))
