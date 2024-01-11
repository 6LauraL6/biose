// MyButton.tsx
import React from 'react';

interface MyButtonProps {
  onClick: () => void;
  buttonText: string;
}

const MyButton: React.FC<MyButtonProps> = ({ onClick, buttonText }) => {
  return (
    <button onClick={onClick}>
      {buttonText}
    </button>
  );
};

export default MyButton;